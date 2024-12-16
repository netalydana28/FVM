#include <bits/stdc++.h>  
using namespace std;
// vector<double> split_double(string line, char delimeter){
//     string my_str = "";
//     vector<double> v;
//     for(int i=0; i<line.size(); i++){
//         if(((line[i] == delimeter) or (i==line.size()-1))){
//             if(line[i+1] == delimeter){
//                 continue;
//             }
//             if(i==line.size()-1){
//                 my_str+=line[i];  
//             }
//             char* temp = new char[my_str.length()+1];
//             strcpy(temp, my_str.c_str());
//             double to_insert = atof(temp);
//             v.push_back(to_insert);
//             my_str = "";
//         }else{
//             // cout<<line[i]<<" here"<<endl;
//             my_str+=line[i];  
//         }
//     }
//     return v;
// }
int main(){
    int n = 81, m = 81, i, j;
    double ro = 1.0, d = 1.0/200.0, dx, dy, dt, f[n][m], newP[n][m], oldP[n][m], newU[n][m], oldU[n][m], newV[n][m], oldV[n][m], e=0.001, max;
    dx = 1.0/(double(n)-1);
    dy = 1.0/(double(m)-1);
    dt = 0.000001;

    for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                oldP[i][j] = 0;
                oldU[i][j] = 0;
                oldV[i][j] = 0;
            newP[i][j] = 0;
                newU[i][j] = 0;
                newV[i][j] = 0;
                f[i][j] = 0;
            }
        }

    for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
        oldU[0][j] = 1;
    }

    for(int k = 0; k<1500000; k++){

        cout<<k<<endl;

	    for(i=0; i<n; i++){
            oldU[i][0] = 0;
            oldU[i][m-1] = 0;
            oldV[i][0] = 0;
            oldV[i][m-1] = 0;
        }
        for(j=0; j<m; j++){
            oldU[0][j] = 0;
            oldU[n-1][j] = newU[n-2][j];
            oldV[0][j] = 0;
            oldV[n-1][j] = newU[n-2][j];
        }
        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            oldU[0][j] = 1;
        }	

        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                newU[i][j] = oldU[i][j] + dt * (d * ((oldU[i+1][j] - 2 * oldU[i][j] + oldU[i-1][j])/(dx*dx)
                    				   + (oldU[i][j+1] - 2 * oldU[i][j] + oldU[i][j-1])/(dy*dy))
                    				    - oldU[i][j] * (oldU[i][j] - oldU[i-1][j])/dx
                    				    - oldV[i][j] * (oldU[i][j] - oldU[i][j-1])/dy);
                newV[i][j] = oldV[i][j] + dt * (d * ((oldV[i+1][j] - 2 * oldV[i][j] + oldV[i-1][j])/(dx*dx)
                    				   + (oldV[i][j+1] - 2 * oldV[i][j] + oldV[i][j-1])/(dy*dy))
                    				    - oldU[i][j] * (oldV[i][j] - oldV[i-1][j])/dx
                    				    - oldV[i][j] * (oldV[i][j] - oldV[i][j-1])/dy);
            }
        }  
	    for(i=0; i<n; i++){
            newU[i][0] = 0;
            newU[i][m-1] = 0;
            newV[i][0] = 0;
            newV[i][m-1] = 0;
        }

        for(j=0; j<m; j++){
            newU[0][j] = 0;
            newU[n-1][j] = newU[n-2][j];
            newV[0][j] = 0;
            newV[n-1][j] = newU[n-2][j];
        }

        for(j=((n-1)/2) - ((n-1)*0.1); j<((n-1)/2) + ((n-1)*0.1) + 1; j++){
            newU[0][j] = 1;
        }

        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                f[i][j] = (ro/dt) * (((newU[i][j] - newU[i-1][j])/dx) + ((newV[i][j] - newV[i][j-1])/dy));
            }
        }

        do{

		    for(i=0; i<n; i++){
            	oldP[i][0] = oldP[i][1];
           	 	oldP[i][m-1] = oldP[i][m-2];
        	}
        	for(j=0; j<m; j++){
           		oldP[0][j] = 0;
            	oldP[n-1][j] = oldP[n-2][j];
       		 }

            for(i=1; i<n-1; i++){
                for(j=1; j<m-1; j++){
                    newP[i][j] = 0.25 * (- (f[i][j] * dx * dx) + oldP[i+1][j] + oldP[i-1][j] + oldP[i][j-1] + oldP[i][j+1]);
                }
            }

            max = 0.0;

            for(i=0; i<n; i++){
                newP[i][0] = newP[i][1];
                newP[i][m-1] = newP[i][m-2];
            }
            for(j=0; j<m; j++){
                newP[0][j] = 0;
                newP[n-1][j] = newP[n-2][j];
            }

            for(i=0; i<n; i++){
                for(j=0; j<m; j++){
                    if(max<fabs(oldP[i][j] - newP[i][j])){
                        max = fabs(oldP[i][j] - newP[i][j]);
                    }
                }
	        }

  	    for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                oldP[i][j] = newP[i][j];
            }
        }
            
        }while(max>e);
       
        for(i=1; i<n-1; i++){
            for(j=1; j<m-1; j++){
                oldU[i][j] = newU[i][j] - ((dt/ro) *(newP[i+1][j] - newP[i][j])/dx);
                oldV[i][j] = newV[i][j] - ((dt/ro) *(newP[i][j+1] - newP[i][j])/dy);           
            }
        }
        if(k == 750000){
            ofstream yFile1("out2M250K.dat");
            yFile1 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<endl;
            for(i=0; i<n; i++){
                for(j=0; j<m; j++){
                    yFile1 << i*dx << " " << j*dy <<" "<<oldP[i][j]<<" "<<oldU[i][j]<<" "<<oldV[i][j]<<endl;
                }
            }
            yFile1.close();
        }
        else if(k == 1000000){
            ofstream yFile2("out2M500K.dat");
            yFile2 << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<endl;
            for(i=0; i<n; i++){
                for(j=0; j<m; j++){
                    yFile2 << i*dx << " " << j*dy <<" "<<oldP[i][j]<<" "<<oldU[i][j]<<" "<<oldV[i][j]<<endl;
                }
            }
            yFile2.close();
        }
    }
    
    ofstream yFile("out3M.dat");
    yFile << "VARIABLES= \"X\", \"Y\", \"P\", \"U\", \"V\"" <<endl;
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            yFile << i*dx << " " << j*dy <<" "<<oldP[i][j]<<" "<<oldU[i][j]<<" "<<oldV[i][j]<<endl;
        }
    }
    yFile.close();
    
    return 0;
}