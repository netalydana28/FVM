#include <bits/stdc++.h> 
#include <sys/time.h> 
#define pi 3.14159265358979323846 
#define w freopen("sol.txt", "w", stdout) 
#define w2 freopen("an.txt", "w", stdout) 
using namespace std; 
const double dx = 1e-3, X = 100.0, alpha = 0.5; 
const int x_size = X / dx + 1; 
double u[x_size], g[x_size]; 
double f(double x) { 
  return x * x; 
} 
double analytical(double x) { 
  return tgamma(2.5 + 1 - 0.5)/tgamma(2.5 + 1) * pow(x, 2.5); 
} 
void preproc() { 
  g[0] = 1; g[1] = -alpha; 
  for (int i = 2; i < x_size; i++) { 
    g[i] = g[i - 1] * (i - 1 - alpha) / i; 
  } 
} 
void solve() { 
 
  struct timeval begin, end; 
  gettimeofday(&begin, 0); 
 
  double maxerr = 0; 
  memset(u, 0, sizeof(u)); 
  u[0] = 0; 
  for (int i = 0; i < x_size - 1; i++) { 
    //printf("%d\n", i); 
    u[i + 1] = f(i* dx)* pow(dx, alpha); 
    for (int k = 1; k <= i + 1; k++) { 
      u[i + 1] -= u[i + 1 - k] * g[k]; 
    } 
    //maxerr = max(maxerr, abs(u[i + 1] - analytical((i + 1) * dx))); 
  } 
 
  gettimeofday(&end, 0); 
  long seconds = end.tv_sec - begin.tv_sec; 
  long microseconds = end.tv_usec - begin.tv_usec; 
  double elapsed = seconds + microseconds * 1e-6; 
  printf("Time measured: %.3f seconds.\n", elapsed); 
  printf("max err in percents: %lf \n", 100 * abs((u[x_size - 1] - analytical((x_size - 1) * dx)) / analytical((x_size - 1) * dx))); 
} 
void show() { 
  w; 
  for (int i = 0; i < x_size; i++) { 
    printf("%lf %lf\n", i * dx, u[i]); 
  } 
} 
void showanalytical() { 
  w2; 
  for (int i = 0; i < x_size; i++) { 
    printf("%lf %lf\n", i * dx, analytical(i * dx)); 
  } 
} 
int main() { 
  preproc(); 
  solve(); 
  show(); 
  showanalytical(); 
}