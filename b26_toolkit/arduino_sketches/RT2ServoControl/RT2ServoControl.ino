#include "Servo.h"

int num_servos = 4;
Servo servos[4];
int servo_ports[] = {5  , 9, 10, 11};

void setup() {
  Serial.begin(9600);
  Serial.setTimeout(1000);
  
  for (int i = 0; i < num_servos; ++i) {
    servos[i].attach(servo_ports[i]);  
  }  
}

void loop() {
  if (Serial.available() > 0) {
    int servo_num = Serial.read();
    while (Serial.available() <= 0) {}
    int pos = Serial.read();
    if (pos) {
      servos[servo_num].write(0);
    }
    else {
      servos[servo_num].write(90);
    }
  }

}
