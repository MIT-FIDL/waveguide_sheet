#include <Wire.h>
#include <Adafruit_ADS1X15.h>

Adafruit_ADS1115 ads_no1, ads_no2, ads_no3, ads_no4;

const int numSensors = 14;
//int16_t adc[numSensors];
long sensorValue[numSensors];
long time_interval = 16; // [millisec]
long time_end = 5000; // [millisec]
long time_start = 0;
long time_base = 0;
long time_now = 0;

void setup(void)
{
  Serial.begin(9600);
// ads1015.setGain(GAIN_TWOTHIRDS);  // 2/3x gain +/- 6.144V  1 bit = 3mV (default)
// ads1115.setGain(GAIN_ONE);     // 1x gain   +/- 4.096V  1 bit = 2mV
// ads1115.setGain(GAIN_FOUR);    // 4x gain   +/- 1.024V  1 bit = 0.5mV
// ads1115.setGain(GAIN_EIGHT);   // 8x gain   +/- 0.512V  1 bit = 0.25mV
// ads1115.setGain(GAIN_SIXTEEN); // 16x gain  +/- 0.256V  1 bit = 0.125mV

  ads_no1.setGain(GAIN_TWO);     // 2x gain   +/- 2.048V  1 bit = 1mV
  ads_no2.setGain(GAIN_TWO);     // 2x gain   +/- 2.048V  1 bit = 1mV
  ads_no3.setGain(GAIN_TWO);     // 2x gain   +/- 2.048V  1 bit = 1mV
  ads_no4.setGain(GAIN_TWO);     // 2x gain   +/- 2.048V  1 bit = 1mV
  ads_no1.begin(0x48);
  ads_no2.begin(0x49);
  ads_no3.begin(0x4A);
  ads_no4.begin(0x4B);
  time_start = micros();
  time_base = micros();
  time_now = micros();
}

void loop(void)
{
//  int16_t adc0, adc1, adc2, adc3;

  for (int i = 0; i < numSensors; i++) {
    sensorValue[i] = get_ads_val(ads_no1, ads_no2, ads_no3, ads_no4, i);
  }
  time_now = micros();

    Serial.print(time_now - time_base);
    for (int i = 0; i < numSensors; i++) {
      Serial.print(",");
      Serial.print(sensorValue[i],9);
    }
    Serial.print(",");
    Serial.println("");
//    delay(2000);
    
}

long get_ads_val(Adafruit_ADS1115 ads_no1,Adafruit_ADS1115 ads_no2,Adafruit_ADS1115 ads_no3,Adafruit_ADS1115 ads_no4,int i)
{
  if (i<=3){
    return ads_no1.readADC_SingleEnded(3-i);
  }
  else if (i<=7){
    return ads_no2.readADC_SingleEnded(i-4);
  }
  else if (i<=11){
    return ads_no3.readADC_SingleEnded(3-(i-8));
  }
  else{
    return ads_no4.readADC_SingleEnded(i-12);
  }
}
