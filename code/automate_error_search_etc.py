from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
import time
import re


# ---------------- helpers ----------------

def min_exposure_time(etc_dict, _f, _snr, _T_eff, _mag):

    _url, _instrument, _filter_set = etc_dict
    _driver = webdriver.Chrome()
    _wait = WebDriverWait(_driver, 20)


    _driver.get(_url)

    

    ###################################################
    # low priority to get the ACS working             #
    # set_filter has to be in a for loop for HST/ACS  #
    # because there are two filter wheels             #
    # instrument and filter would be unique.          #
    # there will be two wheels and two filters for    #
    # each filter (because one will be the clear one) #
    ###################################################

    set_filter(_instrument, _f, _wait)
    # set_signal_to_noise(_snr, _wait)

    # "Time" needed to obtain "SNR" of _snr
    set_exposure(["Time", "SNR"], _snr, _wait)
    
    set_blackbody(_T_eff, _wait)
    
    set_magnitude(_mag, _wait, unit="abmag")
    
    # set_sloan_normalization(_f, _wait)

    set_norm_filter(_f, _filter_set, _wait)

    submit_simulation(_wait)

    request_id, time = extract_time(_wait, _driver)

    # _id_list.append(request_id)

    # if i==0:
    #     keep_it_open = input("wait for me to tell you to close smh")
    
    # exit()
    # keep_it_open = input("wait for me to tell you to close smh")
    _driver.quit()
    return time

def error_lookup(etc_dict, _f, _time, _T_eff, _mags):
    _url, _instrument, _filter_set = etc_dict
    _driver = webdriver.Chrome()
    _wait = WebDriverWait(_driver, 20)


    # for SDSS normalization - i.e., when the magnitudes come from CFHT
    # _norm_filter = "fftype_fNormalizeByFilter.sloan"

    # for HST/WFC3/IR
    # _norm_filter = "fftype_fNormalizeByFilter.wfc3IR"

    _snr_list = []
    for i, T in enumerate(_T_eff):

        _driver.get(_url)

        set_filter(_instrument, _f, _wait)

        set_exposure(["SNR", "Time"], _time, _wait)
        # set_exposure_time(_time, _wait)
        set_blackbody(T, _wait)
        set_magnitude(_mags[i], _wait, unit="abmag")

        set_norm_filter(_f, _filter_set, _wait)
        # set_sloan_normalization(_f, _norm_filter, _wait)



        submit_simulation(_wait)

        request_id, snr = extract_snr(_wait, _driver)

        # _id_list.append(request_id)
        _snr_list.append(snr)

        # if i==0:
        #     keep_it_open = input("wait for me to tell you to close smh")
        
        # exit()
    # keep_it_open = input("wait for me to tell you to close smh")
    _driver.quit()
    return _snr_list




def set_blackbody(temp_k, wait):
    # 1. Select Blackbody radio button
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, "fsorigin_SpectrumBlackBody"
    )))
    radio.click()

    # 2. Wait a moment for JS to activate field
    time.sleep(0.2)

    # 3. Fill temperature field
    temp_field = wait.until(EC.presence_of_element_located((
        By.NAME, "fbbtemp"
    )))

    temp_field.clear()
    temp_field.send_keys(str(int(temp_k)))

def set_filter(instrument, value, wait):

    select_element = wait.until(EC.presence_of_element_located((
        By.NAME, instrument
    )))

    Select(select_element).select_by_value(value)


# def set_filter(selector_names, values, wait):

#     # Single dropdown case
#     if isinstance(selector_names, str):
#         select_element = wait.until(
#             EC.presence_of_element_located(
#                 (By.NAME, selector_names)
#             )
#         )
#         Select(select_element).select_by_value(values)
#         return

#     # Multiple dropdown case
#     elif isinstance(selector_names, (list, tuple)):

#         if not isinstance(values, (list, tuple)):
#             raise ValueError(
#                 "When selector_names is a list, values must also be a list."
#             )

#         if len(selector_names) != len(values):
#             raise ValueError(
#                 "selector_names and values must have the same length."
#             )

#         for selector_name, value in zip(selector_names, values):

#             select_element = wait.until(
#                 EC.presence_of_element_located(
#                     (By.NAME, selector_name)
#                 )
#             )

#             Select(select_element).select_by_value(value)

#         return

#     raise TypeError(
#         "selector_names must be a string or list/tuple of strings."
#     )

def set_exposure(t_or_SNR, val, wait):
    # 1. Select "S/N mode (solve for S/N at given time)"
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, f"simmode_{t_or_SNR[0]}"
    )))
    radio.click()

    # 2. Set exposure time
    time_field = wait.until(EC.presence_of_element_located((
        By.NAME, t_or_SNR[1]#"Time"
    )))
    time_field.clear()
    time_field.send_keys(str(val))



# def set_exposure_time(seconds, wait):
#     # 1. Select "S/N mode (solve for S/N at given time)"
#     radio = wait.until(EC.element_to_be_clickable((
#         By.ID, "simmode_SNR"
#     )))
#     radio.click()

#     # 2. Set exposure time
#     time_field = wait.until(EC.presence_of_element_located((
#         By.NAME, "Time"
#     )))
#     time_field.clear()
#     time_field.send_keys(str(seconds))


# def set_signal_to_noise(snr, wait):
#     # 1. Select "S/N mode (solve for S/N at given time)"
#     radio = wait.until(EC.element_to_be_clickable((
#         By.ID, "simmode_Time"
#     )))
#     radio.click()

#     # 2. Set exposure time
#     time_field = wait.until(EC.presence_of_element_located((
#         By.NAME, "SNR"
#     )))
#     time_field.clear()
#     time_field.send_keys(str(snr))




def set_magnitude(value, wait, unit="abmag"):
    # 1. Select bandpass normalization mode
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, "fftype_fnormalize_bandpass"
    )))
    radio.click()

    # 2. Set flux value
    flux_field = wait.until(EC.presence_of_element_located((
        By.NAME, "rn_flux_bandpass"
    )))
    flux_field.clear()
    flux_field.send_keys(str(value))

    # 3. Select units (AB magnitude etc.)
    unit_dropdown = wait.until(EC.presence_of_element_located((
        By.NAME, "rn_flux_bandpass_units"
    )))

    Select(unit_dropdown).select_by_value(unit)


###
# NEEDS MAJOR IMPROVEMENTS TO APPLY TO HST/WFC3/UVIS AND IR ISOCHRONES
###


def set_norm_filter(filters, filter_set, wait):

    # 1. Select normalization mode
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, f"fftype_fNormalizeByFilter.{filter_set}"
    )))
    radio.click()

    # 2. Map HST filter → Sloan filter or website name 
    mapping = {
        "F336W": "Sloan/U", 
        "F475W": "Sloan/G",
        "F625W": "Sloan/R",
        "F775W": "Sloan/I",
        "F850LP": "Sloan/Z",
        "F105W": "WFC3/IR/F105W", 
        "F110W": "WFC3/IR/F110W", 
        "F125W": "WFC3/IR/F125W", 
        "F140W": "WFC3/IR/F140W",
        "F160W": "WFC3/IR/F160W"
    }

    val = mapping[filters]

    # 3. Set dropdown
    dropdown = wait.until(EC.presence_of_element_located((
        By.NAME, f"filter.{filter_set}" #"filter.wfc3IR" # 
    )))
    Select(dropdown).select_by_value(val)



# def set_sloan_normalization(hst_filter, norm_filter, wait):
#     # 1. Select normalization mode
#     radio = wait.until(EC.element_to_be_clickable((
#         By.ID, norm_filter
#     )))
#     radio.click()

#     # 2. Map HST filter → Sloan filter or website name 
#     mapping = {
#         "F336W": "Sloan/U", 
#         "F475W": "Sloan/G",
#         "F625W": "Sloan/R",
#         "F775W": "Sloan/I",
#         "F850LP": "Sloan/Z",
#         "F105W": "WFC3/IR/F105W", 
#         "F110W": "WFC3/IR/F110W", 
#         "F125W": "WFC3/IR/F125W", 
#         "F140W": "WFC3/IR/F140W",
#         "F160W": "WFC3/IR/F160W"
#     }

#     sloan_value = mapping[hst_filter]

#     # 3. Set dropdown
#     dropdown = wait.until(EC.presence_of_element_located((
#         By.NAME, "filter.sloan" #"filter.wfc3IR" # 
#     )))
#     Select(dropdown).select_by_value(sloan_value)


def submit_simulation(wait):
    btn = wait.until(EC.element_to_be_clickable((
        By.XPATH,
        "//input[@type='button' and @value='Submit Simulation']"
    )))
    btn.click()


def extract_snr(wait, driver):
    # Wait for results page to load (Request ID appears reliably)
    wait.until(EC.presence_of_element_located((
        By.XPATH,
        "//*[contains(text(),'ETC Request ID')]"
    )))

    page_text = driver.page_source

    # Extract ETC Request ID
    request_id = re.search(r"ETC Request ID:\s*([A-Za-z0-9\.\-]+)", page_text)
    request_id = request_id.group(1) if request_id else None

    # Extract SNR value
    snr = re.search(r"SNR\s*=\s*([0-9.]+)", page_text)
    snr = float(snr.group(1)) if snr else None

    return request_id, snr


def extract_time(wait, driver):
    # Wait for results page to load (Request ID appears reliably)
    wait.until(EC.presence_of_element_located((
        By.XPATH,
        "//*[contains(text(),'ETC Request ID')]"
    )))

    page_text = driver.page_source

    # Extract ETC Request ID
    request_id = re.search(r"ETC Request ID:\s*([A-Za-z0-9\.\-]+)", page_text)
    request_id = request_id.group(1) if request_id else None

    # Extract SNR value
    time = re.search(r"Time\s*=\s*([0-9,.]+)", page_text)
    time = time.group(1) if time else None
    time = float(time.replace(",", ""))


    return request_id, time


if __name__ == "__main__":
    

    URL = "https://etc.stsci.edu/etc/input/wfc3uvis/imaging/" #"https://etc.stsci.edu/etc/input/wfc3ir/imaging/" #

    filters = {#{"F105W": "F105W", 
    #             "F110W": "F110W", 
    #             "F125W": "F125W", 
    #             "F140W": "F140W",
    #             "F160W": "F160W"
    #             } 

        "SDSS g": "F475W",
        "SDSS r": "F625W",
        "SDSS i": "F775W",
        "SDSS z": "F850LP"
    }


    exposure_times = {"1hr": 3600} #{
    #     "15min": 900,
    #     "30min": 1800,
    #     "1hr": 3600,
    #     "2hr": 7200,
    #     "5hr": 18000,
    #     "10hr": 36000
    # }

    blackbody_temp = 5000   # K
    ab_mag = 27.0            # AB magnitude
    inst = "wfc3_filter_w" #"irfilt0"#

    driver = webdriver.Chrome()
    wait = WebDriverWait(driver, 20)

    # for HST/WFC3/IR
    # norm_filter = "fftype_fNormalizeByFilter.wfc3IR"
    norm_filter = "fftype_fNormalizeByFilter.sloan"
    # ---------------- main loop ----------------

    results_list = []
    for label, f in filters.items():
        for label, t in exposure_times.items():

            driver.get(URL)

            set_filter(inst, f, wait)
            set_exposure_time(t, wait)
            set_blackbody(blackbody_temp, wait)
            set_normalization(ab_mag, wait, unit="abmag")
            set_sloan_normalization(f, norm_filter, wait)
            
            

            submit_simulation(wait)

            request_id, snr = extract_snr(wait, driver)

            results_list.append({
            	"temp": blackbody_temp,
            	"ab_mag": ab_mag,
                "filter": f,
                # "time_label": label,
                "time_sec": t,
                "snr": snr
            })

            # print(f, label, snr)


    keep_it_open = input("wait for me to tell you to close smh")
    driver.quit()

    print(results_list)
