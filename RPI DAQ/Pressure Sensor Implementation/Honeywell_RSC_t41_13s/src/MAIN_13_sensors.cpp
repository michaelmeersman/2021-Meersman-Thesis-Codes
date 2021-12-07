/* This program uses a Teensy 4.0 to read pressure data from 3 Honeywell RSC sensors */


#include <Arduino.h>
#include "Honeywell_RSC_rework_t41_13s.h"
#include <SPI.h>
#include <Sensor_config_13s.h>


void setup() {
  // Open serial communication
  Serial.begin(115200);

  // Open SPI communication
  SPI.begin();

  // Allowtime to setup SPI
  delay(100);

  // Initialize pressure sensor
  rsc1.init();
  delay(50);
  rsc2.init();
  delay(50);
  rsc3.init();
  delay(50);
  rsc4.init();
  delay(50);
  rsc5.init();
  delay(50);
  rsc6.init();
  delay(50);
  rsc7.init();
  delay(50);
  rsc8.init();
  delay(50);
  rsc9.init();
  delay(50);
  rsc10.init();
  delay(50);
  rsc11.init();
  delay(50);
  rsc12.init();
  delay(50);
  rsc13.init();
  delay(50);

  // Force setup to wait until serial monitor is open (comment if you dont want to use serial monitor!)
  // while (!Serial) {}

  // Print sensor information -- this is not necessary, but is a good sanity check that things are working properly
  // Serial.println();
  // Serial.print("catalog listing:\t");
  // Serial.println(rsc1.catalog_listing());
  // Serial.print("serial number:\t\t");
  // Serial.println(rsc1.serial_number());
  // Serial.print("pressure range:\t\t");
  // Serial.println(rsc1.pressure_range());
  // Serial.print("pressure minimum:\t");
  // Serial.println(rsc1.pressure_minimum());
  // Serial.print("pressure unit:\t\t");
  // Serial.println(rsc1.pressure_unit_name());
  // Serial.print("pressure type:\t\t");
  // Serial.println(rsc1.pressure_type_name());
  // Serial.println();

  delay(50);
}

int multiplier = 1000;
byte buf[26];
short p_vals[13];

void loop() {
  rsc1.get_temperature();
  rsc2.get_temperature(); 
  rsc3.get_temperature();
  rsc4.get_temperature();
  rsc5.get_temperature();
  rsc6.get_temperature();
  rsc7.get_temperature();
  rsc8.get_temperature();
  rsc9.get_temperature();
  rsc10.get_temperature();
  rsc11.get_temperature();
  rsc12.get_temperature();
  rsc13.get_temperature();

  for (int i = 0; i < 10; i++) {

    // rsc1.get_temperature();
    rsc1.select_pressure();
    p_vals[0] = rsc1.read_pressure()*multiplier;

    // rsc2.get_temperature();
    rsc2.select_pressure();
    p_vals[1] = rsc2.read_pressure()*multiplier;

    // rsc3.get_temperature();
    rsc3.select_pressure();
    p_vals[2] = rsc3.read_pressure()*multiplier;

    // rsc4.get_temperature();
    rsc4.select_pressure();
    p_vals[3] = rsc4.read_pressure()*multiplier;

    // rsc5.get_temperature();
    rsc5.select_pressure();
    p_vals[4] = rsc5.read_pressure()*multiplier;

    // rsc6.get_temperature();
    rsc6.select_pressure();
    p_vals[5] = rsc6.read_pressure()*multiplier;

    // rsc7.get_temperature();
    rsc7.select_pressure();
    p_vals[6] = rsc7.read_pressure()*multiplier;

    // rsc8.get_temperature();
    rsc8.select_pressure();
    p_vals[7] = rsc8.read_pressure()*multiplier;

    // rsc9.get_temperature();
    rsc9.select_pressure();
    p_vals[8] = rsc9.read_pressure()*multiplier;

    // rsc10.get_temperature();
    rsc10.select_pressure();
    p_vals[9] = rsc10.read_pressure()*multiplier;

    // rsc11.get_temperature();
    rsc11.select_pressure();
    p_vals[10] = rsc11.read_pressure()*multiplier;

    // rsc12.get_temperature();
    rsc12.select_pressure();
    p_vals[11] = rsc12.read_pressure()*multiplier;

    // rsc13.get_temperature();  
    rsc13.select_pressure();
    p_vals[12] = rsc13.read_pressure()*multiplier;
      
    // for (int i = 0; i<13; i++){
    //   Serial.print(p_vals[i]);
    //   Serial.print('\t');
    // }
    // Serial.println();

    for (int g = 0; g < 13; g++){
      buf[2*g] = p_vals[g] >> 8;
      buf[2*g+1] = p_vals[g] & 0x01ff;
      // Serial.print(p_vals[g]);
      // Serial.print('\t');
    }
    Serial.write(buf,sizeof(buf));
    Serial.println();
  }
}