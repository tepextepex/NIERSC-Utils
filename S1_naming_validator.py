import os.path
import re

"""
Information on Sentinel-1 naming conventions could be found at the ESA web-site:

https://sentinel.esa.int/web/sentinel/user-guides/sentinel-1-sar/naming-conventions
"""


def validate(name, sat=["A", "B"], product=["GRD"], mode=["EW"], res=["H", "M"]):
    """
    Provides a simple check for a file name,
    returning if it is a valid Sentinel-1 archive name.

    Allows to filter products on mission platform, type, sensing mode and resolution class.

    NOTE it does not validate the contents, just a file name!

    @param name: file name to validate "S1...zip"
    @param sat: list of acceptable satellite platforms
    @param product: list of acceptable product types
    @param mode: list of acceptable sensing modes
    @param res: list acceptable resolution classes
    @return: True of False

    """
    satellite = ""
    if "A" in sat:
        satellite += "A"
    if "B" in sat:
        satellite += "B"

    product_type = ""
    if "GRD" in product:
        product_type += "GRD"
    if "RAW" in product:
        product_type += "RAW"
    if "SLC" in product:
        product_type += "SLC"
    if "OCN" in product:
        product_type += "OCN"

    product_type = "".join(set(product_type))  # removes duplicating chars

    resolution_class = ""
    if "H" in res:
        resolution_class += "H"
    if "M" in res:
        resolution_class += "M"
    if "F" in res:
        resolution_class += "F"

    sensor_mode = ""
    if "EW" in mode:
        sensor_mode += "EW"
    if "IW" in mode:
        sensor_mode += "IW"
    if "WV" in mode:
        sensor_mode += "WV"

    sensor_mode = "".join(set(sensor_mode))

    params = {"sat": satellite, "prod": product_type, "res": resolution_class, "mode": sensor_mode}

    s1_pattern = re.compile(
        "S1[%(sat)s]_[%(mode)s]{2}_[%(prod)s]{3}[%(res)s]_1S[SD][HV]_\d{8}T\d{6}_\d{8}T\d{6}_\d{6}_......_.....zip" % params)
    s1_pattern_safe = re.compile(
        "S1[%(sat)s]_[%(mode)s]{2}_[%(prod)s]{3}[%(res)s]_1S[SD][HV]_\d{8}T\d{6}_\d{8}T\d{6}_\d{6}_......_.....SAFE.zip" % params)

    if s1_pattern.match(name) or s1_pattern_safe.match(name):
        return True

    else:
        return False


if __name__ == "__main__":
    print validate("S1B_EW_GRDM_1SDH_20200705T025147_20200705T025228_022330_02A627_3C88.zip")  # True
    print validate(
        "S1B_IW_GRDM_1SDH_20200705T025147_20200705T025228_022330_02A627_3C88.zip")  # False, as IW is rejected by default
    print validate(
        "S1B_EW_GRDF_1SDH_20200705T025147_20200705T025228_022330_02A627_3C88.zip")  # False, as F resolution class is rejected by default
    print validate("S1B_EW_GRDF_1SDH_20200705T025147_20200705T025228_022330_02A627_3C88.zip",
                   res=["F"])  # True, as F resolution class is now included
    print validate("S1B_EW_RAWF_1SDH_20200705T025147_20200705T025228_022330_02A627_3C88.zip",
                   res=["F"])  # False, as RAW product type is not included by default
    print validate("S1B_LOL_KEK_THIS_IS_NOT_A_VALID_NAME_3C88.zip")  # obviously False
