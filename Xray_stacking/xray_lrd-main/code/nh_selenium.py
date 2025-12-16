from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

import numpy as np
from astropy.table import Table

def clean_and_update(ele, val):
    ele.clear()
    ele.send_keys(val)
    
def getnh(ra, dec):
    import subprocess
    
    command = 'nh 2000 ' + str(ra) + ' ' + str(dec)
    print('running ' + command)
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # 检查命令是否成功执行
    if result.returncode == 0:
        output = result.stdout
        # 拆分输出为行
        lines = output.splitlines()
        last_value = None
        # 找到包含 "h1_nh_HI4PI.fits" 的最后一行
        for line in reversed(lines):
            if "h1_nh_HI4PI.fits" in line:
                last_line = line
                # 拆分最后一行并提取最后一个值
                last_value = float(last_line.split()[-1])
                break
        if last_value is None:
            print("Error: No line containing 'h1_nh_HI4PI.fits' found in the output.")
    else:
        print("Error:")
        print(result.stderr)
        last_value = None

    return last_value
def get_nh_milkyway(ra, dec):
    service = Service(ChromeDriverManager().install())
    options = webdriver.ChromeOptions()
    options.add_argument('headless')
    options.add_argument('window-size=1920x1080')
    options.add_argument("disable-gpu")

    driver = webdriver.Chrome(service=service)

    # Navigate to the webpage
    driver.get("https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl")

    # input the coordinate
    clean_and_update(driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/table/tbody/tr[1]/td[2]/input"),"%.5f, %.5f"%(ra, dec))

    driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/p/input[1]").click()

#    clean_and_update(driver.find_element(By.ID, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/table/tbody/tr[6]/td[2]/input"), "0.01")
    clean_and_update(driver.find_element(By.ID, "Radius"), "0.01")

    # Click on the 'Calculate' button
    driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/form/p/input[1]").click()
    driver.implicitly_wait(10)
    # Assuming there is a result field with ID 'result'
    # Retrieve the output value
    output = driver.find_element(By.XPATH, "/html/body/div[2]/div/div[5]/div/div/div[2]/table/tbody/tr/td[2]/center/pre/b[2]")

    textlist = output.text.split(' ')
    print("NH value:", float(textlist[-1]))

    driver.quit()
    return float(textlist[-1])
#    return float(output)

def main():
    # test the code
    tbl_all_lrds = Table.read('/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/narrowlineLRD_cat.fits')
    mwhn_list = []
    for index in range(len(tbl_all_lrds)):
        ra = tbl_all_lrds['RA'][index]
        dec = tbl_all_lrds['Dec'][index]

        print(ra, dec)

        mwhn_list.append(getnh(ra, dec))

    tbl_all_lrds['NH_MW'] = mwhn_list

    tbl_all_lrds.write('/Users/zijianzhang/Astro_Data/LRD_SPEC/Xray_stacking/cstack_stack/narrowlineLRD_cat_mwnh.fits')



if __name__=='__main__':
    main()
