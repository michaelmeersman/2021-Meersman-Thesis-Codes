#!/usr/bin/python3
import serial, string, datetime, time, Adafruit_ADS1x15
import RPi.GPIO as GPIO


xbee = serial.Serial('/dev/ttyS0', 57600)

#Set up signal to NI DAQ
GPIO.setmode(GPIO.BCM)
GPIO.setup(24, GPIO.OUT)


print("waiting for start command...")
xbee.write(str.encode("\r\n"))
pcheck = 0
while True:
	if xbee.in_waiting > 0:
		message = xbee.read(xbee.in_waiting)
		message = message.decode()
		if message == "s": # Or skip decode and check if equal to 73
			print("Trigger on.")
			xbee.write(str.encode("Trigger on."))
			xbee.write(str.encode("\r\n"))
		else:
			time.sleep(0.25)
			break		

		while True:
			# Set GPIO high to indicate start of collection for NI DAQ
			GPIO.output(24,1)
	
			if xbee.read(xbee.in_waiting).decode() == 'e':
				#Set GPIO low to indicate end of collection for NI DAQ
				GPIO.output(24,0)
				print("Trigger off.")
				xbee.write(str.encode("Trigger off."))
				xbee.write(str.encode("\r\n"))
				time.sleep(0.25)
				break
			else:
				continue
	else:
		pcheck += 1
		if pcheck == 10:
			print("Waiting...")
			xbee.write(str.encode("Waiting..."))
			xbee.write(str.encode("\r\n"))
			pcheck = 0
		time.sleep(0.25)
		continue