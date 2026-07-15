from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
import time
import re



# ---------------- helpers ----------------



def error_lookup(_url, _instrument, _f, _time, _T_eff, _mags):

    _driver = webdriver.Chrome()
    _wait = WebDriverWait(_driver, 20)


    _snr_list = []
    for i, T in enumerate(_T_eff):

        _driver.get(_url)

        set_filter(_instrument, _f, _wait)
        set_exposure_time(_time, _wait)
        set_blackbody(T, _wait)
        set_normalization(_mags[i], _wait, unit="abmag")
        set_sloan_normalization(_f, _wait)



        submit_simulation(_wait)

        request_id, snr = extract_results(_wait, _driver)

        # _id_list.append(request_id)
        _snr_list.append(snr)


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

def set_exposure_time(seconds, wait):
    # 1. Select "S/N mode (solve for S/N at given time)"
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, "simmode_SNR"
    )))
    radio.click()

    # 2. Set exposure time
    time_field = wait.until(EC.presence_of_element_located((
        By.NAME, "Time"
    )))
    time_field.clear()
    time_field.send_keys(str(seconds))



def set_normalization(value, wait, unit="abmag"):
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

def set_sloan_normalization(hst_filter, wait):
    # 1. Select Sloan normalization mode
    radio = wait.until(EC.element_to_be_clickable((
        By.ID, "fftype_fNormalizeByFilter.sloan"
    )))
    radio.click()

    # 2. Map HST filter → Sloan filter
    mapping = {
        "F336W": "Sloan/U", 
        "F475W": "Sloan/G",
        "F625W": "Sloan/R",
        "F775W": "Sloan/I",
        "F850LP": "Sloan/Z"
    }

    sloan_value = mapping[hst_filter]

    # 3. Set dropdown
    dropdown = wait.until(EC.presence_of_element_located((
        By.NAME, "filter.sloan"
    )))
    Select(dropdown).select_by_value(sloan_value)


def submit_simulation(wait):
    btn = wait.until(EC.element_to_be_clickable((
        By.XPATH,
        "//input[@type='button' and @value='Submit Simulation']"
    )))
    btn.click()


def extract_results(wait, driver):
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


if __name__ == "__main__":
    

    URL = "https://etc.stsci.edu/etc/input/wfc3uvis/imaging/"

    filters = {
        "SDSS g": "F475W",
        "SDSS r": "F625W",
        "SDSS i": "F775W",
        "SDSS z": "F850LP"
    }


    exposure_times = {
        "15min": 900,
        "30min": 1800,
        "1hr": 3600,
        "2hr": 7200,
        "5hr": 18000,
        "10hr": 36000
    }

    blackbody_temp = 5000   # K
    ab_mag = 27.0            # AB magnitude
    inst = "wfc3_filter_w"

    driver = webdriver.Chrome()
    wait = WebDriverWait(driver, 20)


    # ---------------- main loop ----------------

    results_list = []
    for label, f in filters.items():
        for label, t in exposure_times.items():

            driver.get(URL)

            set_filter(inst, f)
            set_exposure_time(t)
            set_blackbody(blackbody_temp)
            set_normalization(ab_mag, unit="abmag")
            set_sloan_normalization(f)
            
            

            submit_simulation()

            request_id, snr = extract_results()

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
