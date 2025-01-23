import os
import json
import shutil
import logging
import dateparser
import numpy as np
import pandas as pd
import seawater as sw
from copy import deepcopy
from pyrsktools import RSK
import matplotlib.pyplot as plt
from datetime import datetime, timezone
from matplotlib.backends.backend_pdf import PdfPages


def create_file_list(path):
    filetypes = {".tob": "sea&sun", ".rsk": "rbr", ".cnv": "seabird"}
    file_groups = {}
    files = []
    for file in os.listdir(path):
        base_name, extension = os.path.splitext(file)
        if base_name not in file_groups:
            file_groups[base_name] = []
        file_groups[base_name].append(extension)
    for base_name, ext_list in file_groups.items():
        matching_key = next((key for key in ext_list if key.lower() in filetypes), None)
        if matching_key is not None:
            files.append({
                "basename": base_name,
                "extension": matching_key,
                "files": ext_list,
                "type": filetypes[matching_key.lower()],
                "path": os.path.join(path, base_name + matching_key)
            })
    return files


def read_data(file_path, file_type):
    if file_type == "sea&sun":
        profiles = read_sea_and_sun(file_path)
    elif file_type == "rbr":
        profiles = read_rbr(file_path)
    elif file_type == "seabird":
        profiles = read_seabird(file_path)
    else:
        raise ValueError("File type not recognised: {}".format(file_type))
    return profiles


def read_sea_and_sun(file_path):
    column_conversion = {
        "Tur": "Turb",
    }
    skip_rows, columns, units, valid = parse_sea_and_sun(file_path, "Lines :")
    if not valid:
        return False
    df = pd.read_csv(file_path, sep='\s+', header=None, skiprows=skip_rows, names=columns, engine='python',
                     encoding="cp1252")
    columns = list(df.columns)
    for i in range(len(df.columns)):
        if df.columns[i] in column_conversion:
            columns[i] = column_conversion[df.columns[i]]
    df.columns = columns
    df["time"] = parse_time(df)
    df["time"] = df["time"].dt.tz_localize('UTC').astype('int64') // 10 ** 9
    downcast, upcast, air_pressure = extract_single_profile(df)
    profiles = casts_to_profiles(df, downcast, upcast, file_path, "Sea&Sun", air_pressure)
    if len(profiles) == 0:
        return False
    else:
        return profiles


def read_rbr(file_path):
    column_conversion = {
        "timestamp": "time",
        "pressure": "Press",
        "temperature": "Temp",
        "conductivity": "Cond",
        "chlorophyll": "Chl_A",
        "dissolved_o2_saturation": "sat",
        "dissolved_o2_concentration": "DO_mg",
    }
    air_pressure = False
    with RSK(file_path) as rsk:
        rsk.readdata()
        df = pd.DataFrame(rsk.data["timestamp"], columns=["time"])
        for column in rsk.channelNames:
            if column in column_conversion:
                df[column_conversion[column]] = rsk.data[column]
        df["time"] = df["time"].dt.tz_localize('UTC').astype(int) // 10 ** 3

        rsk.computeprofiles()
        downcast = rsk.getprofilesindices(direction="down")
        upcast = rsk.getprofilesindices(direction="up")
        air_idx = np.setdiff1d(np.arange(len(df)), np.array(flatten([downcast, upcast])))
        if len(air_idx) > 0:
            pressure = np.array(df["Press"])
            air_pressure = np.nanmean(pressure[air_idx])

    profiles = casts_to_profiles(df, downcast, upcast, file_path, "RBR", air_pressure)
    if len(profiles) == 0:
        return False
    else:
        return profiles


def read_seabird(file_path):
    column_conversion = {
        "Pressure, Digiquartz [db]": "Press",
        "Pressure, Strain Gauge [db]": "Press",
        "Temperature [ITS-90, deg C]": "Temp",
        "Conductivity [uS/cm]": "Cond",
        "Turbidity, WET Labs ECO [NTU]": "Turb",
        "Oxygen, SBE 43 [% saturation]": "sat",
        "Oxygen, SBE 43 [mg/l]": "DO_mg",
        "Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]": "Flur"
    }
    data_conversion = {
        "Conductivity [uS/cm]": divide_by_1000,
    }
    skip_rows, columns, valid, time = parse_seabird(file_path, "*END*")
    if not valid:
        return False
    df = pd.read_csv(file_path, sep='\s+', header=None, skiprows=skip_rows, names=columns, engine='python',
                     encoding="cp1252")
    df["time"] = df["Time, Elapsed [seconds]"] + time

    columns = list(df.columns)
    for i in range(len(df.columns)):
        if df.columns[i] in data_conversion:
            df[df.columns[i]] = data_conversion[df.columns[i]](df[df.columns[i]])
        if df.columns[i] in column_conversion:
            columns[i] = column_conversion[df.columns[i]]
    df.columns = columns

    downcast, upcast, air_pressure = extract_single_profile(df)
    profiles = casts_to_profiles(df, downcast, upcast, file_path, "Seabird", air_pressure)
    if len(profiles) == 0:
        return False
    else:
        return profiles


def casts_to_profiles(df, downcast, upcast, file_path, file_type, air_pressure):
    profiles = []
    for index, cast in enumerate(downcast):
        if index == 0:
            name = os.path.splitext(os.path.basename(file_path))[0]
        else:
            name = os.path.splitext(os.path.basename(file_path))[0] + "_{}".format(index)
        bottom = cast[-1] - cast[0]
        df_profile = df.iloc[cast[0]: upcast[index][-1] + 1]
        data = {
            "name": name,
            "file": os.path.basename(file_path),
            "folder": os.path.dirname(file_path),
            "type": file_type,
            "bottom_index": bottom,
            "data": df_profile
        }
        if air_pressure:
            data["air_pressure"] = air_pressure
        profiles.append(data)
    return profiles


def parse_sea_and_sun(input_file_path, string):
    valid = True
    with open(input_file_path, encoding="latin1", errors='ignore') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if string in lines[i]:
            break
    columns = lines[i + 2].replace(";", "").split()
    columns.pop(0)
    columns = rename_duplicates(columns)
    units = lines[i + 3].replace(";", "").replace("[", "").replace("]", "").replace("°", "deg").split()
    skip_rows = i + 5
    n = 0
    while len(lines[i + 5].split()) - 1 > len(columns):
        columns.append(n)
        n = n + 1
    if len(lines) <= skip_rows + 1 or len(columns) < 5:
        valid = False
    return skip_rows, columns, units, valid


def parse_seabird(input_file_path, string):
    valid = True
    with open(input_file_path, encoding="latin1", errors='ignore') as f:
        lines = f.readlines()
    columns = []
    for i in range(len(lines)):
        if "# start_time" in lines[i]:
            time = dateparser.parse(lines[i].split("= ")[1].split(" [")[0])
            time = time.replace(tzinfo=timezone.utc).timestamp()
        if "# name" in lines[i]:
            columns.append(lines[i].split(":")[1].lstrip().replace("\n", ""))
        if string in lines[i]:
            break
    skip_rows = i + 1
    if len(lines) <= skip_rows + 1 or len(columns) < 5:
        valid = False
    return skip_rows, columns, valid, time


def process_profiles(profiles, folder, template):
    os.makedirs(folder, exist_ok=True)
    pdf_pages = PdfPages(os.path.join(folder, 'profiles_requiring_metadata.pdf'))
    plots_per_page = 6
    total_plots = len(profiles)

    for i in range(0, total_plots, plots_per_page):
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(8.5, 11))
        axes = axes.flatten()

        for j in range(plots_per_page):
            plot_index = i + j
            if plot_index < total_plots:
                profile = profiles[plot_index]
                time = datetime.fromtimestamp(np.array(profile["data"]["time"])[0]).strftime('%H:%M %d %B %Y')
                x = np.array(profile["data"]["Temp"])
                y = np.array(profile["data"]["Press"]) * -1
                axes[j].plot(x, y, color="lightgrey")
                x[profile["bottom_index"] + 1:-1] = np.nan
                axes[j].plot(x, y, color="red", label=profile["type"])
                axes[j].set_title(profile["name"] + ".meta\n" + time)
                axes[j].set_xlabel("Temperature (°C)")
                axes[j].set_ylabel("Pressure (dbar)")
                axes[j].legend()
            else:
                axes[j].axis('off')
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        pdf_pages.savefig(fig)
        plt.close(fig)
    pdf_pages.close()

    with open(template, 'r') as f:
        metadata = json.load(f)
    for profile in profiles:
        m = metadata.copy()
        time = datetime.fromtimestamp(np.array(profile["data"]["time"])[0])
        m["filename"] = profile["file"]
        m["campaign"]["Device"] = profile["type"]
        m["campaign"]["Date of measurement"] = time.strftime('%Y-%m-%d')
        m["profile"]["Profile name"] = profile["name"]
        m["profile"]["Time of measurement (local)"] = time.strftime('%H:%M')
        with open(os.path.join(folder, profile["name"] + ".meta"), 'w') as f:
            json.dump(m, f, indent=4)


def parse_time(df):
    if "IntD" in df.columns and "IntT" in df.columns:
        return pd.to_datetime(df["IntD"] + " " + df["IntT"], format="%d.%m.%Y %H:%M:%S.%f", dayfirst=True)
    elif "IntDT" in df.columns and "IntDT1" in df.columns:
        return pd.to_datetime(df["IntDT"] + " " + df["IntDT1"], format="%d.%m.%Y %H:%M:%S.%f", dayfirst=True)
    elif "IntT" in df.columns and "IntT1" in df.columns:
        return pd.to_datetime(df["IntT"] + " " + df["IntT1"], format="%d.%m.%Y %H:%M:%S.%f", dayfirst=True)


def extract_single_profile(df, rolling=3, diff=0.01, var="Cond", pressure="Press", max_pressure_cut=3.0):
    start_index = 0
    pressure_arr = np.array(df[pressure])
    var_arr = np.array(df[var])
    air_pressure = False
    try:
        max_start = np.min(pressure_arr) + max_pressure_cut
        df = pd.DataFrame(var_arr)
        df_mean = df.rolling(rolling, center=True).mean().bfill().ffill()
        df_diff = df_mean.diff()
        outliers = np.array(np.abs(df_diff) > diff).flatten()
        start_index = np.where(outliers[:np.argmax(pressure_arr > max_start)])[0][-1]
    except:
        logging.warning("Failed to locate start of profile")
    bottom_index = np.argmax(pressure_arr)
    downcast = [list(range(start_index, bottom_index + 1))]
    upcast = [list(range(bottom_index + 1, len(df)))]
    if start_index > 0:
        air_pressure = np.nanmean(pressure_arr[:start_index])
    return downcast, upcast, air_pressure


def divide_by_1000(arr):
    return arr / 1000


def flatten(nested_list):
    flattened = []
    for item in nested_list:
        if isinstance(item, list):
            flattened.extend(flatten(item))
        else:
            flattened.append(item)
    return flattened


def rename_duplicates(arr):
    out = []
    d = {}
    for i in arr:
        d.setdefault(i, -1)
        d[i] += 1
        if d[i] >= 1:
            out.append('%s%d' % (i, d[i]))
        else:
            out.append(i)
    return out


def ch1903_to_latlng(x, y):
    x_aux = (x - 600000) / 1000000
    y_aux = (y - 200000) / 1000000
    lat = 16.9023892 + 3.238272 * y_aux - 0.270978 * x_aux ** 2 - 0.002528 * y_aux ** 2 - 0.0447 * x_aux ** 2 * y_aux - 0.014 * y_aux ** 3
    lng = 2.6779094 + 4.728982 * x_aux + 0.791484 * x_aux * y_aux + 0.1306 * x_aux * y_aux ** 2 - 0.0436 * x_aux ** 3
    lat = (lat * 100) / 36
    lng = (lng * 100) / 36
    return lat, lng


def first_centered_differences(x, y, fill=False):
    if x.size != y.size:
        raise ValueError("Vectors do not have the same size")
    dy = np.full(x.size, np.nan)
    iif = np.where((np.isfinite(x)) & (np.isfinite(y)))[0]
    if iif.size == 0:
        return dy
    x0 = x[iif]
    y0 = y[iif]
    dy0 = np.full(x0.size, np.nan)
    # calculates differences
    dy0[0] = (y0[1] - y0[0]) / (x0[1] - x0[0])
    dy0[-1] = (y0[-1] - y0[-2]) / (x0[-1] - x0[-2])
    dy0[1:-1] = (y0[2:] - y0[0:-2]) / (x0[2:] - x0[0:-2])

    dy[iif] = dy0

    if fill:
        dy[0:iif[0]] = dy[iif[0]]
        dy[iif[-1] + 1:] = dy[iif[-1]]
    return dy


def json_converter(qa):
    for keys in qa.keys():
        try:
            if qa[keys]["simple"]["bounds"][0] == "-inf":
                qa[keys]["simple"]["bounds"][0] = -np.inf
            if qa[keys]["simple"]["bounds"][1] == "inf":
                qa[keys]["simple"]["bounds"][1] = np.inf
        except:
            pass
    try:
        if qa["time"]["simple"]["bounds"][1] == "now":
            qa["time"]["simple"]["bounds"][1] = datetime.now().timestamp()
        return qa
    except:
        return qa


def copy_files(file_path, folder, extensions):
    directory, filename = os.path.split(file_path)
    filename_no_ext = os.path.splitext(filename)[0]
    files_in_directory = os.listdir(directory)
    for file in files_in_directory:
        name, extension = os.path.splitext(file)
        if filename_no_ext in name:
            src_path = os.path.join(directory, file)
            dest_path = os.path.join(folder, file)
            if extension.lower() in extensions:
                shutil.copy2(src_path, dest_path)


def copy_variables(variables_dict):
    var_dict = dict()
    for var in variables_dict:
        var_dict[var] = variables_dict[var][:]
    nc_copy = deepcopy(var_dict)
    return nc_copy

def default_salinity_temperature(temperature):
    return 1.8626 - 0.052908 * temperature + 0.00093057 * temperature ** 2 - 6.78e-6 * temperature ** 3

def salinity(Temp, Cond, y_cond, temperature_func= default_salinity_temperature):
    ft = temperature_func(Temp)
    cond20 = ft * Cond * 1000
    salin = y_cond * cond20
    return salin

def density(temperature, salinity):
    rho = 1e3 * (
                0.9998395 + 6.7914e-5 * temperature - 9.0894e-6 * temperature ** 2 + 1.0171e-7 * temperature ** 3 -
                1.2846e-9 * temperature ** 4 + 1.1592e-11 * temperature ** 5 - 5.0125e-14 * temperature ** 6 + (
                    8.181e-4 - 3.85e-6 * temperature + 4.96e-8 * temperature ** 2) * salinity)
    return rho

def oxygen_saturation(T, S, altitude=372., lat=46.2, units="mgl"):
    # calculates oxygen saturation in mg/l according to Garcia-Benson
    # to be coherent with Hannah
    if units != "mgl" and units != "mll":
        units = "mgl"
    mgL_mlL = 1.42905
    mmHg_mb = 0.750061683
    mmHg_inHg = 25.3970886
    standard_pressure_sea_level = 29.92126
    standard_temperature_sea_level = 15 + 273.15
    gravitational_acceleration = gr = sw.g(lat)
    air_molar_mass = 0.0289644
    universal_gas_constant = 8.31447
    baro = (1. / mmHg_mb) * mmHg_inHg * standard_pressure_sea_level * np.exp(
        (-gravitational_acceleration * air_molar_mass * altitude) / (
                    universal_gas_constant * standard_temperature_sea_level))
    u = 10 ** (8.10765 - 1750.286 / (235 + T))
    press_corr = (baro * mmHg_mb - u) / (760 - u)

    Ts = np.log((298.15 - T) / (273.15 + T))
    lnC = 2.00907 + 3.22014 * Ts + 4.0501 * Ts ** 2 + 4.94457 * Ts ** 3 + -0.256847 * Ts ** 4 + 3.88767 * Ts ** 5 - S * (
                0.00624523 + 0.00737614 * Ts + 0.010341 * Ts ** 2 + 0.00817083 * Ts ** 3) - 4.88682e-07 * S ** 2
    O2sat = np.exp(lnC)
    if units == "mll":
        O2sat = O2sat * press_corr
    elif units == "mgl":
        O2sat = O2sat * mgL_mlL * press_corr

    return O2sat

def potential_temperature_sw(T, S, p, p_ref):
    """
    Calculates potential temperature as per UNESCO 1983 report.
    Parameters
    ----------
    s(p) : array_like
        salinity [psu (PSS-78)]
    t(p) : array_like
        temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    pr : array_like
        reference pressure [db], default = 0
    Returns
    -------
    pt : array_like
        potential temperature relative to PR [℃ (ITS-90)]
    """
    return sw.ptmp(s=S,t=T,p=p,pr=p_ref)

def position_in_array(arr, value):
    for i in range(len(arr)):
        if value < arr[i]:
            return i
    return len(arr)