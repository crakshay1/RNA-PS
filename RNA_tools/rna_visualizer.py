import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

def represent(path_to_file):
    driver = webdriver.Chrome()
    driver.get("https://molstar.org/viewer")

    # Getting rid of a button
    driver.find_element(By.CSS_SELECTOR,"button[title='Load a structure from the provided source and create its representation.']").click()
    # Finding the file input area
    driver.find_element(By.CSS_SELECTOR,"button[title='Load one or more files and optionally create default visuals']").click()
    file_input = driver.find_element(By.CSS_SELECTOR, "input[type='file']")
    file_input.send_keys(os.path.abspath(path_to_file))
    driver.find_element(By.CSS_SELECTOR, "button[class='msp-btn msp-btn-block msp-btn-commit msp-btn-commit-on']").click()