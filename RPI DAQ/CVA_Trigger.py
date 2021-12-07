import serial, string, datetime, time, Adafruit_ADS1x15, math
import RPi.GPIO as GPIO

#Set up signal to NI DAQ
GPIO.setmode(GPIO.BCM)
GPIO.setup(24, GPIO.OUT)

GPIO.output(24,1)

i = 0
while True:
    i+=1
    