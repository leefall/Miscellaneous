#!/usr/bin/env python
from selenium import webdriver
from selenium.webdriver.common.by import By
import time
import random
import os
StartTime=(time.ctime())
def enable_download_headless(browser,download_dir):
	browser.command_executor._commands["send_command"] = ("POST", '/session/$sessionId/chromium/send_command')
	params = {'cmd':'Page.setDownloadBehavior', 'params': {'behavior': 'allow', 'downloadPath': download_dir}}
	browser.execute("send_command", params)



def Crawl(sFileID):
	
	print(sFileID)
	
	try:
		os.mkdir("/mnt/alpha/leefall2/TCGA_LGG/Metadata/"+sFileID)
	except:
		pass
	
	
	options = webdriver.ChromeOptions()
	options.add_argument('--no-sandbox')
	options.add_argument('--disable-gpu')
	prefs={"download.default_directory":"/mnt/alpha/leefall2/TCGA_LGG/Metadata/"+sFileID,"download.prompt_for_download":False,
	"download.directory_upgrade": True,
	"safebrowsing_for_trusted_sources_enabled": False,
	"safebrowsing.enabled": False
	}
	
	options.add_experimental_option("prefs",prefs)
	options.add_argument('--verbose')
	options.add_argument('--disable-software-rasterizer')
	options.add_argument('--headless')
	
		
	
	sDriver="/storage/home/leefall2/tools/chromedriver/chromedriver"
	#sWebsite="https://portal.gdc.cancer.gov/legacy-archive/files/db9c1922-bf47-4dff-ab62-2964aeaef587"
	sWebsite="https://portal.gdc.cancer.gov/legacy-archive/files/"+sFileID
	driver = webdriver.Chrome(sDriver, options=options)
	driver.set_window_size(1900, 1200)
	
	enable_download_headless(driver, "/mnt/alpha/leefall2/TCGA_LGG/Metadata/"+sFileID)
	
	
	driver.get(sWebsite)
	nRand=random.randrange(30,50)
	driver.implicitly_wait(time_to_wait=nRand)
	try:
		downloadzip=driver.find_element(By.CLASS_NAME,'btn-primary')
		downloadzip.click()
		nRand=random.randrange(10,15)
		time.sleep(nRand)
		
		
		
		downloadzip=driver.find_elements(By.CLASS_NAME,'btn-primary')
		
		for i in range(2,len(downloadzip)):
			downloadzip[i].click()
			nRand=random.randrange(6,15)
			time.sleep(nRand)
	
	except:
		print("Error~!~!~!~!~!~!~!~!")
		driver.save_screenshot(sFileID+'.png')
	
	
	driver.close()
	

#fp=open("/mnt/alpha/leefall2/TCGA_OV/CrongReview_NoC239_TCGAaliquotID_table.txt")
fp=open("/mnt/alpha/leefall2/TCGA_LGG/ForReview_mc3_NonNA_gdc_manifest_20211116_094540.txt")
fp.readline()
#fp.readline()
for sLine in fp.readlines():
	t=sLine.split("\t")
	Crawl(t[0])
	#break

#Crawl("db9c1922-bf47-4dff-ab62-2964aeaef587")




print(StartTime)
print(time.ctime())
