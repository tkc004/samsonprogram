import os, errno

def TKmkdir(dl_path):
    if not os.path.exists(dl_path):
        print "path doesn't exist. trying to make"
        os.makedirs(dl_path)