//unsigned long next = 0;

void setup() {
  Serial.begin(115200);
}

void loop() {
  int sensorValue = analogRead(A0);
  Serial.println(sensorValue);
  delay(1);
  //  if (micros() >= next) {
    //next = micros() + 1000;           // one sample every 1 ms = 1000 Hz
    //Serial.println(analogRead(A0));   // MyoWare output (RAW or ENV) wired to A0
  //}
}