// Define all of the pin allocations for the LP
const int LP_SWEEP = 35;
const int LP_POL = 30;
const int LP_EN = 29;
const int LP_I_POS = 14;
const int LP_I_NEG = 12;
const int LP_MOTOR = 23;

// Physical constants (in SI units)
const float k_B = 1.38064852E-23; // Boltzmann's constant
const float Ap = 1.0; // Probe surface area
const float m_i = 1; // Mass of the ions
const float q_e = 1.60217662E-19; // Elementary charge

// Current sensor calibration constants
double IPOS_OFFSET = 0.0; // Default value
double INEG_OFFSET = 0.0; // Default value
const double CURRENT_SCALE_POS = double(2.0); // 2 for positive branch
const double CURRENT_SCALE_NEG = double(1.0/13.47); // 1/13.47 for negative branch

// Other constants
const int NUM_SAMPLES = 250; // Nummber of samples to average over for a single data point
const int delay_sample = 0;  // Delay between sampling data during the sweep
const int fit_size = 20; // number of data points to use for the linear fit

// Function to calibrate sensor offset
void CalibrateLP(){

  const int num = 500;

  int sample_p;
  int sample_n;

  IPOS_OFFSET = 0;
  INEG_OFFSET = 0;

  digitalWrite(LP_EN, LOW); // Enable the probes
  digitalWrite(LP_POL, HIGH); // Forward polarity
  analogWrite(LP_SWEEP, 0);  // Write the probe voltage to zero
  delay(100);

  for (int i = 0; i < num; i++){
    sample_p = analogRead(LP_I_POS);
    sample_n = analogRead(LP_I_NEG);
    IPOS_OFFSET += double(sample_p);
    INEG_OFFSET += double(sample_n);
    
  }

  // Average the results
  IPOS_OFFSET = IPOS_OFFSET/double(num);
  INEG_OFFSET = INEG_OFFSET/double(num);

  // Convert to voltage CHANGE FOR TEENSY
  IPOS_OFFSET *= 5.0/1023.0;
  INEG_OFFSET *= 5.0/1023.0;

  digitalWrite(LP_EN, HIGH); // Disable the probes
  
}

// Function to sample a single data point
// Returns a voltage in V
double AverageSample(int num, bool positive){
  double res_p = 0;
  double res_n = 0;
  double res;
  int sample_p;
  int sample_n;

  for (int i = 0; i < num; i++){
    sample_p = analogRead(LP_I_POS);
    sample_n = analogRead(LP_I_NEG);
    res_p += double(sample_p);
    res_n += double(sample_n);
    
  }

  // Average the results
  res_p = res_p/double(num);
  res_n = res_n/double(num);

  // Convert to voltage CHANGE FOR TEENSY
  res_p *= 5.0/1023.0;
  res_n *= 5.0/1023.0;

  // Apply the offset
  res_p -= IPOS_OFFSET;
  res_n -= INEG_OFFSET;

  // Convert to single voltage
  res = res_p + res_n;

  // Apply the scale factor
  if (positive == true){
    res *= CURRENT_SCALE_POS;
  }
  else
  {
    res *= CURRENT_SCALE_NEG;
  }
  
  return (res);

}

// Function for performing the sweep and returning the data
void sweep(double result[511]){
  
  digitalWrite(LP_EN, LOW); // Enable the probes
  digitalWrite(LP_POL, LOW); // Reverse polarity
  
  // Perform the negative sweep first
  for (int i = -255; i < 0; i++){
    analogWrite(LP_SWEEP, -i);  // Write the probe voltage
    delay(delay_sample);    // Wait before sampling
    result[i+255] = -AverageSample(NUM_SAMPLES, false); // Read the sensor
  }
  
  digitalWrite(LP_POL, HIGH); // Forward polarity
  
  // Perform the positve sweep second
  for (int i = 0; i <= 255; i++){
    analogWrite(LP_SWEEP, i); // Write the probe voltage
    delay(delay_sample);    // Wait before sampling
    result[i+255] = AverageSample(NUM_SAMPLES, true); // Read the sensor
  }
  
  // Set the output of the probes to zero
  analogWrite(LP_SWEEP, 0);
  
  digitalWrite(LP_EN, HIGH); // Disable the probes
}

// Function for initialising the LP pins
void LP_Setup(){
  
  pinMode(LP_POL, OUTPUT);
  pinMode(LP_EN, OUTPUT);
  pinMode(LP_SWEEP, OUTPUT);
  
}

// Function to analyse the data obtained from the probe
void LP_analyse(double data[511], double V[511], double properties[2]){

  // Initialise some variables
  double ne; // electron density
  double Te; // electron temperature
  double fit_is[2]; // ion saturation region fit array
  double fit_origin[2]; // origin region fit array
  double I_is; // electron saturation current
  double slope_origin; // slope about the origin

  // 1. Linear fit in ion saturation region
  double V_is[fit_size];
  double data_is[fit_size];
  ArrayExtract(V,V_is,0,fit_size);
  ArrayExtract(data,data_is,0,fit_size);
  simpLinReg(V_is, data_is, fit_is, fit_size);
  I_is = fit_is[1];

  // 2. Linear fit about the origin
  double V_origin[2*fit_size+1];
  double data_origin[2*fit_size+1];
  ArrayExtract(V,V_origin,256-fit_size,256+fit_size+1);
  ArrayExtract(data,data_origin,256-fit_size,256+fit_size+1);
  simpLinReg(V_origin,data_origin, fit_origin, 2*fit_size+1);
  slope_origin = fit_origin[0];  

  // 3. Compute the electron temperature
  Te = (q_e*I_is)/(2*k_B*slope_origin);
  
  // 4. Compute the electron density
  ne = I_is/(0.61*q_e*Ap*sqrt(k_B*Te/m_i));
  
  // Return the plasma properties
  properties[0] = ne;
  properties[1] = Te;
}

// Linear regression function
void simpLinReg(double x[], double y[], double lrCoef[2], int n){
// pass x and y arrays (pointers), lrCoef pointer, and n. The lrCoef array is comprised of the slope=lrCoef[0] and intercept=lrCoef[1]. n is length of the x and y arrays.
// http://en.wikipedia.org/wiki/Simple_linear_regression

// initialize variables
double xbar=0;
double ybar=0;
double xybar=0;
double xsqbar=0;

// calculations required for linear regression
for (int i=0; i<n; i++){
xbar=xbar+x[i];
ybar=ybar+y[i];
xybar=xybar+x[i]*y[i];
xsqbar=xsqbar+x[i]*x[i];
}
xbar=xbar/n;
ybar=ybar/n;
xybar=xybar/n;
xsqbar=xsqbar/n;

// simple linear regression algorithm
lrCoef[0]=1; // (xybar-xbar*ybar)/(xsqbar-xbar*xbar);
lrCoef[1]=ybar-lrCoef[0]*xbar;
}

// Function to generate the voltage array
void getVoltage(double V[], double Vmax, int n){

  int i;
  for (i = -n; i < (n+1); i++) {
    V[i+n] = Vmax*double(i) / n;

  }
}

// Function to extract one array out of another array
void ArrayExtract(double original[], double output[], int a, int b){

  for(int i = 0 ; i < (b-a) ; ++i){
    output[i] = original[i+a];

  }
}
