#include "Servo.h"

Servo servo1;

int pos = 0;

void setup() {
  Serial.begin(9600);
  Serial.setTimeout(1000);
  
  servo1.attach(9);

  
}

void loop() {
  if (Serial.available() > 0) {
    int servoNum = Serial.read();
    while (Serial.available() <= 0) {}
    int pos = Serial.read();
    if (pos) {
      servo1.write(0);
    }
    else {
      servo1.write(90);
    }
  }

}
