import os
import os.path
from S1L1Tools import S1L1Tools
from S1_naming_validator import validate
import glob  # TODO: really, you are able to find all the *.xmls without glob!


def render_directory(source_dir, out_dir=None):
    if out_dir is None:
        out_dir = os.path.join(source_dir, "renders")  # if doesn't exist, will be created automatically
    source_dir = os.path.join(source_dir, "")  # adds a trailing "/", if omitted

    if os.path.isdir(source_dir):  # checks also if exists
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            print "%s directory created" % out_dir

        for file_name in os.listdir(source_dir):
            full_path = os.path.join(source_dir, file_name)
            if os.path.isfile(full_path) and validate(file_name, mode=["EW", "IW"]):
                s = S1L1Tools(full_path)
                print "Rendering scene %s into %s" % (file_name, out_dir)
                s.render(out_dir, ['HH'], preserve_source_file_name=True)
                del s

        xmls = glob.glob(os.path.join(out_dir, "*.tif.aux.xml"))
        for xml in xmls:
            os.remove(xml)


if __name__ == "__main__":
    source_dir = "/home/tepex/NIERSC/IEPI/20200709/source/"  # directory with S1 zip-archives. May contain other files - they won't be processed
    render_directory(source_dir)
