#!/usr/bin/python3
import serial, string, datetime, time, Adafruit_ADS1x15
import RPi.GPIO as GPIO


#Set up signal to NI DAQ
GPIO.setmode(GPIO.BCM)
GPIO.setup(24, GPIO.OUT)

while True:
	trigger = int(input('0 for trigger off, 1 for trigger on:'))
	if trigger == 0:
		GPIO.output(24,0)
		print("trigger is off...")
	if trigger == 1:
		GPIO.output(24,1)
		print("trigger is on...")
	else:
		continue