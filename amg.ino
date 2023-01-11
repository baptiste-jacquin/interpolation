#include <Wire.h>
#include <Adafruit_AMG88xx.h>

Adafruit_AMG88xx amg;
float pixels[AMG88xx_PIXEL_ARRAY_SIZE];

void setup() {
    Serial.begin(2000000);
    bool status;
    status = amg.begin();
    delay(100);
}

void loop() { 
    amg.readPixels(pixels);
    for(int i=0; i<AMG88xx_PIXEL_ARRAY_SIZE; i++){
      Serial.print(pixels[i]);
      Serial.print(" ");
    }
    Serial.println();
    delay(100);
}