#include <bits/stdc++.h>  
using namespace std;

// #define n 51 
// #define m 51 

int main(){
    double u = 1, v = 0, ro = 1, d=1.0/500.0;
    int i, j, iter = 0; 
    int n = 81, m = 81;
    double dx, dy, dt, f[n][m], newP[n][m], oldP[n][m], oldU[n][m], newU[n][m], oldV[n][m], l1[n][m], l2[n][m], u_star[n][m], v_star[n][m], newV[n][m], eps, max; 
    dx = 1/(double(n)-1); 
    dy = 1/(double(m)-1);
    dt = 0.000088;
    eps = 0.00001;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            oldP[i][j] = 0; 
            oldU[i][j] = 0; 
            oldV[i][j] = 0;
            newU[i][j] = 0;
            newV[i][j] = 0; 
            l1[i][j] = 0;
            l2[i][j] = 0;
            u_star[i][j] = 0;
            v_star[i][j] = 0;
            f[i][j] = 0;
        }
    } 
    for(int k=0; k<5000; k++){
        cout << k << endl;
        for(j=0; j<m; j++){
            oldU[0][j] = 0;
            oldV[0][j] = 0;
            oldU[n-1][j] = oldU[n-2][j]; 
            oldV[n-1][j] = oldV[n-2][j];
        }
        for(i=0; i<n; i++){
            oldU[i][m-1] = 0;
            oldU[i][0] = 0;
            oldV[i][m-1] = 0;
            oldV[i][0] = 0;
        }
        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            oldU[0][j] = 1;
        } 
        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                l1[i][j] = d * (((oldU[i+1][j] - 2*oldU[i][j] + oldU[i-1][j])/(dx*dx)) 
                               +((oldU[i][j+1]-2*oldU[i][j]+oldU[i][j-1])/(dy*dy))) 
                               - oldU[i][j]*((oldU[i][j] - oldU[i-1][j])/dx) 
                               - oldV[i][j]*((oldU[i][j] - oldU[i][j-1])/dy);
                l2[i][j] = d * (((oldV[i+1][j] - 2*oldV[i][j] + oldV[i-1][j])/(dx*dx)) 
                              + ((oldV[i][j+1]-2*oldV[i][j]+oldV[i][j-1])/(dy*dy)))
                              - oldU[i][j]*((oldV[i][j] - oldV[i-1][j])/dx) 
                              - oldV[i][j]*((oldV[i][j] - oldV[i][j-1])/dy);
            }
        }
        for(j=0; j<m; j++){
            l1[0][j] = 0;
            l2[0][j] = 0;
            l1[n-1][j] = l1[n-2][j]; 
            l2[n-1][j] = l2[n-2][j];
        }
        for(i=0; i<n; i++){
            l1[i][m-1] = 0;
            l1[i][0] = 0;
            l2[i][m-1] = 0;
            l2[i][0] = 0;
        }
        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            l1[0][j] = 1;
        } 
        // for(j = 0; j<m; j++){
        //     u_star[n-1][j] = u_star[n-2][j]; 
        //     v_star[n-1][j] = v_star[n-2][j];
        // }
        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                u_star[i][j] = l1[i][j]*dt + oldU[i][j];
                v_star[i][j] = l2[i][j]*dt + oldV[i][j];
            }
        }
        for(j=0; j<m; j++){
            u_star[0][j] = 0;
            v_star[0][j] = 0;
            u_star[n-1][j] = u_star[n-2][j]; 
            v_star[n-1][j] = v_star[n-2][j];
        }
        for(i=0; i<n; i++){
            u_star[i][m-1] = 0;
            u_star[i][0] = 0;
            v_star[i][m-1] = 0;
            v_star[i][0] = 0;
        }
        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            u_star[0][j] = 1;
        } 

        for(int r = 0; r < 100; r++){
            for(i=1; i<n-1; i++){
                for(j=1; j<m-1; j++){
                    f[i][j] = (ro/dt) * (((u_star[i][j] - u_star[i-1][j])/dx) + ((v_star[i][j] - v_star[i][j-1])/dy));
                }
            }
            for(j = 0; j<m; j++){
                oldP[0][j] = 0;
                oldP[n-1][j] = oldP[n-2][j]; 
            }
            for(i=0; i<n; i++){
                oldP[i][0] = oldP[i][1];
                oldP[i][m-1] = oldP[i][m-2];
            }
            for(i=1; i<n-1; i++){
                for(j=1; j<m-1; j++){
                    newP[i][j] = 0.25*(-(f[i][j] * dx * dx) + oldP[i+1][j] + oldP[i-1][j] + oldP[i][j+1] + oldP[i][j-1]);
                }
            }
            for(i=0; i<n; i++){
                for(j=0; j<m; j++){
                    oldP[i][j] = newP[i][j];
                }
            }
            for(j = 0; j<m; j++){
                oldP[0][j] = 0;
                oldP[n-1][j] = oldP[n-2][j]; 
            }
            for(i=0; i<n; i++){
                oldP[i][0] = oldP[i][1];
                oldP[i][m-1] = oldP[i][m-2];
            }
            for(j = 0; j<m; j++){
                newP[0][j] = 0;
                newP[n-1][j] = newP[n-2][j];
            }
            for(i=0; i<n; i++){
                newP[i][0] = newP[i][1];
                newP[i][m-1] = newP[i][m-2];
            }
        }
        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                newU[i][j] = u_star[i][j] - ((dt/ro)*((newP[i+1][j] - newP[i][j])/dx));
                newV[i][j] = v_star[i][j] - ((dt/ro)*((newP[i][j+1] - newP[i][j])/dy));
            }
        }
        for(j=0; j<m; j++){
            newU[0][j] = 0;
            newV[0][j] = 0;
            newU[n-1][j] = newU[n-2][j]; 
            newV[n-1][j] = newV[n-2][j];
        }
        for(i=0; i<n; i++){
            newU[i][m-1] = 0;
            newU[i][0] = 0;
            newV[i][m-1] = 0;
            newV[i][0] = 0;
        }
        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            newU[0][j] = 1;
        } 
        for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                oldP[i][j] = newP[i][j];
                oldU[i][j] = newU[i][j];
                oldV[i][j] = newV[i][j];
            }

        }
       
    }
    ofstream yFile("out6.dat");
    yFile << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<endl;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            yFile << i*dx << " " << j*dy <<" "<<newP[i][j]<<" "<<newU[i][j]<<" "<<newV[i][j]<<endl;
        }
    }
    yFile.close();
    return 0;
}