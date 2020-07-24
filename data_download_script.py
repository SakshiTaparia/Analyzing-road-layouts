
#THE SCRIPT for automating the data retrieval


import requests
import time
#the user should ensure that (maxlon-minlon)/dx and (maxlat-minlat)/dy are integers

dx = 0.01 #float(input())
dy = 0.01 #float(input())

minlon = 72.755 #float(input())
minlat = 18.895 #float(input())
maxlon = 73.095 #float(input())
maxlat = 19.295 #float(input())

cur_minlon = minlon
cur_minlat = minlat
cur_maxlon = minlon + dx
cur_maxlat = minlat + dy

it = 136
length = round(maxlon*100 - minlon*100)
print(length)
while(cur_maxlat<=maxlat+0.000000001):
  cur_minlon = minlon
  cur_maxlon = minlon + dx
  while(cur_maxlon<=maxlon+0.000000001):
    time.sleep(10)
    #print(f"{cur_minlon:.4f}" + ',' + f"{cur_minlat:.4f}" + ',' + f"{cur_maxlon:.4f}" + ',' + f"{cur_maxlat:.4f}")
    print("attempting ",end="")
    print(it+1)
    print("To resume:",end="")
    print("it = " + str((int(it/length)*length)),end=" ")
    
    #Change the constant value as per the original min_lat of every district
    print("minlat = " + str(float(int(it/length))*0.01 + 18.855))
    my_referer = "https://www.openstreetmap.org/export"
    
    res = requests.get('https://api.openstreetmap.org/api/0.6/map?bbox='+ f"{cur_minlon:.4f}" + ',' + f"{cur_minlat:.4f}" + ',' + f"{cur_maxlon:.4f}" + ',' + f"{cur_maxlat:.4f}", headers = {'referer': my_referer})
    
    #res = requests.get('http://overpass.openstreetmap.ru/cgi/xapi_meta?*[bbox='+ f"{cur_minlon:.4f}" + ',' + f"{cur_minlat:.4f}" + ',' + f"{cur_maxlon:.4f}" + ',' + f"{cur_maxlat:.4f}" + ']', headers = {'referer': my_referer})
    if res.status_code != 200:
      continue
    found = False
    for chunk in res.iter_content(1000):
      c = str(chunk)
      if (c.find("You have downloaded too much data")!=-1):
        found=True
      break;
      
    if(found):
      continue
    it = it + 1
    playFile = open('data' + str(it) + '.txt', 'wb')
    for chunk in res.iter_content(100000):
            playFile.write(chunk)
    cur_maxlon = cur_maxlon + dx
    cur_minlon = cur_minlon + dx
  cur_maxlat = cur_maxlat + dy
  cur_minlat = cur_minlat + dy
  
  
  
  
'''
from google.colab import files

i = 1
while i<=it:
  files.download('data' + str(i) + '.txt')
'''
