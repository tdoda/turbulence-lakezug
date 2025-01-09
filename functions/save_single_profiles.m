function msg=save_single_profiles(BINNED, SLOW, FAST, fname)

save(fname,'BINNED','SLOW','FAST')
msg='ok';