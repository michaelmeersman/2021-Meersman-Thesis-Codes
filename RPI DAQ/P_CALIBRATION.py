#!/usr/bin/python3
import serial, string, datetime, time, Adafruit_ADS1x15, struct, adafruit_gps
import RPi.GPIO as GPIO
import numpy as np

# Wrapper function for pyserial readline function
class ReadLine:
	def __init__(self,s):
		self.buf = bytearray()
		self.s = s
		# print(s)

	def readline(self):
		# MY Readline Function
		# self.buf = bytearray()
		while True:
			i = self.buf.find(b"\n")
			if i >= 0:
				r = self.buf[:i-1]
				self.buf = self.buf[i+1:]
				if len(r) == 26:
					return r
			data = self.s.read(60) # if I always read twice the length of a frame of data, then I guarantee I will have at lease 1 full frame
			i = data.find(b"\n")
			if i >= 0: # check if I found a new line character in the bytes I read
				q = self.buf + data[:i-1] # by choosing k-1 I eliminate both \r\n bytes
				self.buf[0:] = data[i+1:]
				if len(q) == 26: # if my current byte array, 'r' is the correct length print it out
					return q
			else:
				self.buf.extend(data)

xbee = serial.Serial('/dev/ttyS0',57600)
ser1 = serial.Serial('/dev/serial/by-path/platform-fd500000.pcie-pci-0000:01:00.0-usb-0:1.3:1.0',115200)
ser2 = serial.Serial('/dev/serial/by-path/platform-fd500000.pcie-pci-0000:01:00.0-usb-0:1.4:1.0',115200)
ser3 = serial.Serial('/dev/serial/by-path/platform-fd500000.pcie-pci-0000:01:00.0-usb-0:1.1:1.0',115200)
ser4 = serial.Serial('/dev/serial/by-path/platform-fd500000.pcie-pci-0000:01:00.0-usb-0:1.2:1.0',115200)

# Create object for readline wrapper function above
r1 = ReadLine(ser1)
r2 = ReadLine(ser2)
r3 = ReadLine(ser3)
r4 = ReadLine(ser4)

# Create file string for calibration data
current_date = input('Input current date (mm_dd_yy): ')
flight = input('Input flight number (Flight_#): ')
location = "/home/pi/Documents/data_" + current_date + "_" + flight + "/"
filename = "CALIBRATION_" + current_date + "_" + flight + ".txt"
calfilestring = location + filename
calfile = open(calfilestring,'w')

print("Beginning calibration sequence...")
xbee.write(str.encode("----------------------------------------------------------------------------------\r\n"))
xbee.write(str.encode("Beginning calibration sequence...\r\n"))

cal_sample = int(input('Input number of calibration samples to collect:'))
cal = [[0 for x in range(52)] for y in range(cal_sample)]
p = 1
for k in range(cal_sample):
	output1 = r1.readline()
	output2 = r2.readline()
	output3 = r3.readline()
	output4 = r4.readline()
	output = output1+output2+output3+output4
	print(p)			
	output = struct.unpack('>hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh', output)
	cal[k][0:52] = output[0:52]
	p += 1

cal_array = np.array(cal)
CAL = np.mean(cal_array, axis=0)
print(CAL)
calfloat = [0 for x in range(52)]
calstr = ''
for i in range(52):
	calfloat[i] = float(CAL[i])
	calstr += '{0:.4f}'.format(CAL[i]) + '\t'

calfile.write(calstr)
calfile.close()
print("Calibration complete!")
xbee.write(str.encode("Calibration complete!\r\n"))