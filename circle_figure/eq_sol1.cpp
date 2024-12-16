//Circle with triangle mesh
//Working with input data
#include <bits/stdc++.h>  
using namespace std;
const int node_num=1407;
map<int, vector<double> > nodes;
const int center_num=2656;
map<int, vector<double> > centers;
int face_num;
map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_cell_to_node;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<int, vector<int> > node_to_cell;
map<vector<int> , double> sn;
vector<int> boundary_face;
vector<int> nodes_boundary;
double P_nodes[node_num];
vector<double> P_faces;
double P_centers_old[center_num];
double P_centers_new[center_num], max_diff;
const double e=0.000001;
double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(int i; i<line.size(); i++){
        if(((line[i] == delimeter) or (i==line.size()-1))){
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            double to_insert = atof(temp);
            v.push_back(to_insert);
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}
vector<int> split_int(string line, char delimeter){
    string my_str = "";
    vector<int> v;
    for(int i; i<line.size(); i++){
        // cout<<my_str<<" "<<endl;
        if((line[i] == delimeter) or (i==line.size()-1)){
            if(line[i+1] == delimeter){
                continue;
            }
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            if(my_str[0]== 'C'){
                my_str = "";
                continue;
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            int to_insert = atof(temp);
            v.push_back(to_insert);
            // cout<<"Vector: "<< v[v.size()-1]<<endl;
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}
void nodes_coor(){
    string myText;
    ifstream NodesFile("sq_coor.txt");
    vector<double> v; 
    int cnt = 1;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        // cout<<v[0]<<" "<<v[1]<<endl;
        v.erase(v.end()-1);
        nodes[cnt] = v;
        v.clear();
        cnt++;
    }
    NodesFile.close();
}
void cell_to_node(){
    string myText;
    ifstream NodesFile("sq_links.txt");
    vector<int> v; 
    vector<int> v1; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(int i=1; i<4; i++){
            v1.push_back(cell);
            v1.push_back(i);
            l_cell_to_node[v1] = v[i-1];
            if(node_to_cell.find(v[i-1]) != node_to_cell.end()){
               node_to_cell[v[i-1]].push_back(cell);
            }else{
                vector<int> temp;
                temp.push_back(cell);
                node_to_cell[v[i-1]] = temp;
                temp.clear();
            }
            v1.clear();
        }
        v.clear();
    }
    NodesFile.close();
}
void center_coordinates_generator(){
    vector<int> v;
    vector<double> v1;
    for(int i = 1; i < center_num + 1; i++){
        v.push_back(i);
        v.push_back(1);
        register int node1 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(2);
        register int node2 = l_cell_to_node[v];
        v.erase(v.end()-1);
        v.push_back(3);
        register int node3 = l_cell_to_node[v];
        v1.push_back((nodes[node1][0] + nodes[node2][0] + nodes[node3][0]) / 3);
        v1.push_back((nodes[node1][1] + nodes[node2][1] + nodes[node3][1]) / 3);
        centers[i]=v1;
        v.clear();
        v1.clear();
    }
}
void face_links(){
    map< set<int>, vector<int> > nodes_to_cell;
    vector<int> key;
    vector<int> key1;
    register int cnt = 1;
    vector<int> trg_nodes;
    set<int> nodes_ind;
    vector< vector<int> > c_t_f(center_num);
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<4; j++){
            key.push_back(j);
            register int trg_node = l_cell_to_node[key];
            trg_nodes.push_back(trg_node);
            key.erase(key.end()-1);
        }
        key.clear();
        register int m;
        for(int j=1; j<4; j++){
            if(j==3){
                m = 0;
            }else{
                m = j;
            }
            nodes_ind.insert(trg_nodes[j-1]); 
            nodes_ind.insert(trg_nodes[m]);
            if(nodes_to_cell.find(nodes_ind) == nodes_to_cell.end()){
                key1.push_back(i);
                nodes_to_cell[nodes_ind] = key1; 
                key1.clear(); 

            }else{
                nodes_to_cell[nodes_ind].push_back(i);
            }
            nodes_ind.clear();           
        }
        trg_nodes.clear();
    }
    face_num = nodes_to_cell.size();
    map< set<int>, vector<int> >::iterator m_it=nodes_to_cell.begin();
    register int f_cnt = 1;
    while(m_it!=nodes_to_cell.end()){
        key1 = m_it->second; 
        set <int>::iterator s_it = m_it->first.begin();
        register int s_cnt = 1;
        while(s_it != m_it->first.end()){
            key.push_back(f_cnt); 
            key.push_back(s_cnt);
            l_face_to_node[key] = *s_it;
            s_cnt++;
            s_it++;
            key.clear();
        }
        key.push_back(f_cnt);
        for(int i=1; i<key1.size()+1; i++){
            c_t_f[key1[i-1]-1].push_back(f_cnt);
            key.push_back(i);
            l_face_to_cell[key] = key1[i-1]; 
            key.erase(key.end()-1);
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                nodes_boundary.push_back(l_face_to_node[key]);
                key.push_back(i+1);
                // nodes_boundary.push_back(l_face_to_node[key]);
                l_face_to_cell[key] = key1[i-1];
                key.erase(key.end()-1);
            }
        }
        f_cnt++;
        m_it++;
        key1.clear();
        key.clear();
    }

    // for(int i=0; i<center_num; i++){
    //     for(int j = 0; j<c_t_f[i].size(); j++){
    //         cout<<c_t_f[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    // cout<<c_t_f[0].size()<<endl;
    key.push_back(boundary_face[1]);
    key.push_back(1);
    key1.push_back(boundary_face[2]);
    key1.push_back(0);
    if(l_face_to_node[key] != l_face_to_node[key1]){
        int temp = boundary_face[0];
        boundary_face[0] = boundary_face[1];
        boundary_face[1] = temp;
    }
    key.clear();
    key1.clear();
    for(int i=1; i<center_num+1; i++){
        key.push_back(i);
        for(int j=1; j<c_t_f[i-1].size()+1; j++){
            key.push_back(j);
            l_cell_to_face[key] = c_t_f[i-1][j-1];
            // cout<<l_cell_to_face[key]<<" "<<endl;
            key.erase(key.end()-1);
        }
        key.clear();
    }
}
bool is_boundary_face(int face_ind){
    vector<int> key1;
    key1.push_back(face_ind);
    key1.push_back(1);
    vector<int> key2; 
    key2.push_back(face_ind);
    key2.push_back(2);
    if(l_face_to_cell[key1] == l_face_to_cell[key2]){
        return true;
    }
    return false;
}
vector<double> face_center(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    vector<double> ans; 
    ans.push_back((n1[0] + n2[0])/2);
    ans.push_back((n1[1] + n2[1])/2);
    return ans;
}
double find_max(double  v[center_num]){
    double max = v[1];
    for(int i=1; i<center_num; i++){
        if(v[i]>max){
            max = v[i];
        }
    }
    return max;
}
vector<double> normal(int face_ind){
    vector <double> ans; 
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    ans.push_back(n1[1] - n2[1]); 
    ans.push_back(n2[0] - n1[0]);  
    double len = pow(pow(ans[0], 2) + pow(ans[1], 2), 0.5); 
    for(int i=0; i<ans.size(); i++){
        ans[i] /= len;
    }
    return ans; 
}
double area(int face_ind){
    vector <int> key; 
    key.push_back(face_ind); 
    key.push_back(1);
    vector<double> n1 = nodes[l_face_to_node[key]];
    key.erase(key.end()-1);
    key.push_back(2);
    vector<double> n2 = nodes[l_face_to_node[key]];
    return pow((pow((n1[0] - n2[0]), 2) + pow((n1[1] - n2[1]), 2)) , 0.5);
}
double cell_to_vertex_interpolation(int node_ind){ 
    vector<double> node_coor = nodes[node_ind]; 
    vector<double> center_coor;  
    int cell_ind; 
    double nom = 0; 
    double denom = 0; 
    for(int i=0; i<node_to_cell[node_ind].size(); i++){ 
        cell_ind = node_to_cell[node_ind][i]; 
        center_coor = centers[cell_ind]; 
        double d = pow(pow((center_coor[0] - node_coor[0]), 2) + pow((center_coor[1] - node_coor[1]), 2), 0.5); 
        double P = P_centers_old[cell_ind - 1]; 
        nom += P/d; 
        denom += 1/d; 
    } 
    return nom/denom; 
}
void update_nodes(){
    // register int bound1 = 2;//bottom
    // register int bound2 = 21;
    for(int i=1; i<node_num+1; i++){
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()){
                continue;
        }else{
            P_nodes[i-1] = cell_to_vertex_interpolation(i);
        }
    }
    vector<int> v; 
    for(int i=0; i<boundary_face.size(); i++){
        // P_faces[121-1] = 8;
        int b_f = boundary_face[i];
        v.push_back(b_f); 
        v.push_back(1);
        int n1 = l_face_to_node[v];
        double x1 = nodes[n1][0];
        double y1 = nodes[n1][1];
        // P_nodes[l_face_to_node[v]-1] = 8; 
        v.erase(v.end()-1);
        v.push_back(2);
        int n2 = l_face_to_node[v];
        double x2 = nodes[n2][0];
        double y2 = nodes[n2][1];
        if((x1 == 0) or (x1 == 0.25) or (y1 == 0.5)){
            P_nodes[n1-1] = (sin(x1+2*y1))+(exp(2*x1+3*y1));
            P_nodes[n2-1] = (sin(x2+2*y2))+(exp(2*x2+3*y2));
        }
        // if((x2 == 0) or (x2 == 0.25) or (y2 == 0.5)){'
            
        // }else{
        //     P_nodes[n2-1] = (2*cos(x2+2*y2))+(3*exp(2*x2+3*y2));
        //     nodes_boundary.push_back(n2);
        // }
        P_faces[b_f-1] = (P_nodes[n1-1] + P_nodes[n2-1])/2;
        // P_nodes[l_face_to_node[v]-1 ] = 8; 
        v.clear();
    }
}
void fvm(){
    int iter=1;
    do{
        cout<<"iteration: "<<iter<<"\n";
        iter++;
        for(int i=1; i<center_num+1; i++){
            double denom = 0; 
            double nom = 0;
            vector<int> key;
            key.push_back(i);
            vector<double> center_coor = centers[i];
            cout<<"here"<<endl;
            for(int j=1; j<4; j++){
                key.push_back(j);
                register int face_ind = l_cell_to_face[key];
                vector<int> key1;
                key1.push_back(face_ind);
                double deltx, delty, nx, ny, valid_denom;
                vector<int> v; 
                if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
                    vector<double> face_coor = face_center(face_ind);
                    deltx = face_coor[0] - center_coor[0];
                    delty = face_coor[1] - center_coor[1];
                    key1.push_back(1);
                    nx = sn[key1];
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    ny = sn[key1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    valid_denom = 0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    nom += P_faces[face_ind-1] * area(face_ind) * valid_denom;
                    denom += area(face_ind) * valid_denom;
                    v.push_back(face_ind); 
                    v.push_back(1);
                    int n[2];
                    n[0] = l_face_to_node[v];
                    double x1 = nodes[n[0]][0];
                    double y1 = nodes[n[0]][1];
                    v.erase(v.end()-1);
                    v.push_back(2);
                    n[1] = l_face_to_node[v];
                    v.clear();
                    double x2 = nodes[n[1]][0];
                    double y2 = nodes[n[1]][1];
                    if(fabs(y1 - 0) == 0 and fabs(y2 - 0) == 0){
                        for(int k = 0; k<2; k++){
                            for(int i=0; i<node_to_cell[n[k]].size(); i++){
                                v.push_back(node_to_cell[n[k]][i]);
                                // cout<<"node# "<<n[k]<<" center: "<<node_to_cell[n[k]][i]<<"\n";
                                center_coor = centers[i];
                                // center_coor[0] = centers[node_to_cell[n[k]][i]][0];
                                // center_coor[1] = centers[node_to_cell[n[k]][i]][1];
                                double nom = 0, denom = 0, valid_denom;
                                for(int j=1; j<4; j++){
                                    v.push_back(j);
                                    // cout<<"key "<<v[0]<<" "<<v[1]<<"\n";
                                    register int face_ind = l_cell_to_face[v];
                                    v.erase(v.end()-1);
                                    // cout<<node_to_cell[n[k]][i]<<" "<<j<<" "<<face_ind<<" index"<<"\n";
                                    key.push_back(face_ind);
                                    key.push_back(1);
                                    int fn1 = l_face_to_node[key];
                                    key.erase(key.end()-1);
                                    key.push_back(2);
                                    int fn2 = l_face_to_node[key];
                                    key.clear();
                                    // cout<<"face nodes "<<fn1<<" "<<fn2<<"\n";
                                    if(fn2 == n[k]){
                                        fn2 = fn1; 
                                    }else if ((fn2 != n[k]) and (fn1 != n[k])){
                                        // cout<<"Faceno "<<face_ind<<"\n";
                                        continue;
                                    }
                                    // cout<<"Faceyes "<<face_ind<<"\n";
                                    vector<double> face_coor = face_center(face_ind);
                                    double nx = face_coor[1] - center_coor[1];
                                    double ny = center_coor[0] - face_coor[0];
                                    nx /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                                    ny /= pow(pow(nx, 2) + pow(ny, 2), 0.5);
                                    double area = pow(pow(face_coor[0] - center_coor[0], 2) + pow(face_coor[1] - center_coor[1], 2), 0.5);
                                    double dx = nodes[fn2][0] - nodes[n[k]][0];
                                    double dy = nodes[fn2][1] - nodes[n[k]][1];
                                    if(((dx>0) and (nx<0)) or ((dx<0) and (nx>0))){
                                        nx*=(-1);
                                    }
                                    if(((dy>0) and (ny<0)) or ((dy<0) and (ny>0))){
                                        ny*=(-1);
                                    }
                                    valid_denom=0;
                                    if(dx != 0){
                                        valid_denom += (nx/dx);
                                    }
                                    if(dy != 0){
                                        valid_denom += (ny/dy);
                                    }

                                    nom += P_nodes[fn2-1] * area * valid_denom;
                                    denom += area * valid_denom;
                                    // cout<<" end "<<v.size()<<"\n";
                                    // cout<<"Area "<<area<<" valid denom "<<valid_denom<<"\n";
                                    // cout<<"FACE " <<j<<" "<<nom<<" "<<denom<<"\n";
                                }
                                P_nodes[n[k]-1] = nom/denom;
                                // cout<<" final "<<nom<<" "<<denom<<"\n";
                                v.clear();

                            }

                        }
                    }

                    key1.clear();
                }else{

                    key1.push_back(1);
                    nx = sn[key1];
                    register int out_center_ind = l_face_to_cell[key1]; 
                    key1.erase(key1.end()-1);
                    key1.push_back(2);
                    vector<double> c2 = centers[out_center_ind];
                    ny = sn[key1];

                    if(out_center_ind == i){
                        out_center_ind = l_face_to_cell[key1];
                        c2 = centers[out_center_ind];
                    }
                    deltx = c2[0] - center_coor[0];
                    delty = c2[1] - center_coor[1];
                    if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
                        nx*=(-1);
                    }
                    if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
                        ny*=(-1);
                    }
                    //***
                    valid_denom = 0;
                    if(deltx != 0){
                        valid_denom += (nx/deltx);
                    }
                    if(delty != 0){
                        valid_denom += (ny/delty);
                    }
                    nom += P_centers_old[out_center_ind-1] * area(face_ind) * valid_denom; 
                    denom += area(face_ind) * valid_denom; 
                    key1.clear();
                }
                key.erase(key.end()-1);
            }
            nom += sin(center_coor[0] + 2 * center_coor[1]) - 16 * exp(2 * center_coor[0] + 3 * center_coor[1]);
            denom -= 3;
            P_centers_new[i-1] = nom/denom; 
            key.clear();
        }
        cout<<"here "<<endl;
        max_diff = 0.0; 
        for(int i=0; i<center_num; i++){
            if(max_diff<fabs(P_centers_old[i] - P_centers_new[i])){
                max_diff = fabs(P_centers_old[i] - P_centers_new[i]);
            }
        }
        for(int i=0; i<center_num; i++){
            P_centers_old[i] = P_centers_new[i];
        }
        // vector<int> v; 
        // for(int i=0; i<boundary_face.size(); i++){
        //     int b_f = boundary_face[i];
        //     v.push_back(b_f); 
        //     v.push_back(1);
        //     int n1 = l_face_to_node[v];
        //     double x1 = nodes[n1][0];
        //     double y1 = nodes[n1][1];
        //     v.erase(v.end()-1);
        //     v.push_back(2);
        //     int n2 = l_face_to_node[v];
        //     double x2 = nodes[n2][0];
        //     double y2 = nodes[n2][1];
        //     if((x1 == 0) or (x1 == 0.25) or (y1 == 0.5)){
        //         P_nodes[n1-1] = (sin(x1+2*y1))+(exp(2*x1+3*y1));
        //         P_nodes[n2-1] = (sin(x2+2*y2))+(exp(2*x2+3*y2));
        //     }
        //     P_faces[b_f-1] = (P_nodes[n1-1] + P_nodes[n2-1])/2;
        //     v.clear();
        // }
        update_nodes();
    }while(max_diff>e);
}
// void fvm(){
//     for(int i=1; i<center_num+1; i++){
//         vector<double> gDiff;
//         vector<int> key;
//         key.push_back(i);
//         double ac=0;
//         double bc=0;
//         vector<double> center_coor = centers[i];
//         for(int j=1; j<4; j++){
//             key.push_back(j);
//             register int face_ind = l_cell_to_face[key];
//             vector<int> key1;
//             key1.push_back(face_ind);
//             double deltx;
//             double delty;
//             double nx; 
//             double ny;
//             double valid_denom;
//             if(find(boundary_face.begin(), boundary_face.end(), face_ind) != boundary_face.end()){   
//                 vector<double> face_coor = face_center(face_ind);
//                 deltx = face_coor[0] - center_coor[0];
//                 delty = face_coor[1] - center_coor[1];
//                 key1.push_back(1);
//                 nx = sn[key1];
//                 key1.erase(key1.end()-1);
//                 key1.push_back(2);
//                 ny = sn[key1];
//                 if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
//                     nx*=(-1);
//                 }
//                 if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
//                     ny*=(-1);
//                 }
//                 valid_denom = 0;
//                 if(deltx != 0){
//                     valid_denom += (nx/deltx);
//                 }
//                 if(delty != 0){
//                     valid_denom += (ny/delty);
//                 }
//                 bc += (area(face_ind))*P_faces[face_ind-1]*valid_denom;
//                 ac += (area(face_ind))*valid_denom;
//                 gDiff.push_back((-1)*(area(face_ind))*valid_denom);
//                 // if((deltx != 0) and (delty != 0)){
//                 //     bc += (area(face_ind))*P_faces[face_ind-1]*((nx/deltx) +(ny/delty));
//                 //     ac += (area(face_ind))*((nx/deltx) + (ny/delty));
//                 //     gDiff.push_back((-1)*(area(face_ind))*((nx/deltx) + (ny/delty)));
//                 // }
//                 // else if (deltx == 0){
//                 //     bc += (area(face_ind))*P_faces[face_ind-1]*(ny/delty);
//                 //     ac += (area(face_ind))*(ny/delty);
//                 //     gDiff.push_back((-1)*(area(face_ind))*(ny/delty));
//                 // }else if (delty == 0){
//                 //     bc += (area(face_ind))*P_faces[face_ind-1]*(nx/deltx);
//                 //     ac += (area(face_ind))*(nx/deltx);
//                 //     gDiff.push_back((-1)*(area(face_ind))*(nx/deltx));
//                 // }
//                 key1.clear();
//             }else{
//                 key1.push_back(1);
//                 nx = sn[key1];
//                 register int out_center_ind = l_face_to_cell[key1]; 
//                 key1.erase(key1.end()-1);
//                 key1.push_back(2);
//                 vector<double> c2 = centers[out_center_ind];
//                 ny = sn[key1];

//                 if(out_center_ind == i){
//                     out_center_ind = l_face_to_cell[key1];
//                     c2 = centers[out_center_ind];
//                 }
//                 deltx = c2[0] - center_coor[0];
//                 delty = c2[1] - center_coor[1];
//                 if(((deltx>0) and (nx<0)) or ((deltx<0) and (nx>0))){
//                     nx*=(-1);
//                 }
//                 if(((delty>0) and (ny<0)) or ((delty<0) and (ny>0))){
//                     ny*=(-1);
//                 }
//                 //***
//                 valid_denom = 0;
//                 if(deltx != 0){
//                     valid_denom += (nx/deltx);
//                 }
//                 if(delty != 0){
//                     valid_denom += (ny/delty);
//                 }
//                 gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*valid_denom;
//                 ac += (area(face_ind))*valid_denom;
//                 gDiff.push_back((-1)*(area(face_ind))*valid_denom);

//                 key1.clear();

//             }
//             key.erase(key.end()-1);
//         }
//         gauss_seidel_coeff[i-1][i-1] = ac;
//         gauss_seidel_b[i-1] = bc;
//         ac = 0;
//         bc = 0; 
//         key.clear();
//         gDiff.clear();
//     }
//     gauss_seidel();
// }
void set_boundary_conditions(){
    // register int bound1 = 2;//bottom
    // register int bound2 = 21;

    vector<int> v; 
    for(int i=0; i<boundary_face.size(); i++){
        // P_faces[121-1] = 8;
        int b_f = boundary_face[i];
        v.push_back(b_f); 
        v.push_back(1);
        int n1 = l_face_to_node[v];
        double x1 = nodes[n1][0];
        double y1 = nodes[n1][1];
        // P_nodes[l_face_to_node[v]-1] = 8; 
        v.erase(v.end()-1);
        v.push_back(2);
        int n2 = l_face_to_node[v];
        double x2 = nodes[n2][0];
        double y2 = nodes[n2][1];
        if((x1 == 0) or (x1 == 0.25) or (y1 == 0.5)){
            P_nodes[n1-1] = (sin(x1+2*y1))+(exp(2*x1+3*y1));
            P_nodes[n2-1] = (sin(x2+2*y2))+(exp(2*x2+3*y2));
        }else{
            P_nodes[n1-1] = (2*cos(x1+2*y1))+(3*exp(2*x1+3*y1));
            P_nodes[n2-1] = (2*cos(x1+2*y1))+(3*exp(2*x1+3*y1));
        }

        // if((x2 == 0) or (x2 == 0.25) or (y2 == 0.5)){'
            // nodes_boundary.push_back(n1);
            // nodes_boundary.push_back(n2);
        // }else{
        //     P_nodes[n2-1] = (2*cos(x2+2*y2))+(3*exp(2*x2+3*y2));
        //     nodes_boundary.push_back(n2);
        // }
        P_faces[b_f-1] = (P_nodes[n1-1] + P_nodes[n2-1])/2;
        // P_nodes[l_face_to_node[v]-1 ] = 8; 
        v.clear();
    }
}

int main(){
    cout<<"0"<<endl;
    nodes_coor();
    cout<<"1"<<endl;
    cell_to_node();
    cout<<"2"<<endl;
    center_coordinates_generator();
    cout<<"3"<<endl;
    face_links();
    cout<<"4"<<endl;
    vector<double> v1;
    vector<int> v1_;
    for(int i=1; i<face_num+1; i++){
        v1 = normal(i);
        v1_.push_back(i);
        v1_.push_back(1);
        sn[v1_] = v1[0];
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        sn[v1_] = v1[1];
        v1_.clear();
        v1.clear();
    }
    vector<int> m1;
    vector<int> m2;
    for(int i=0; i<boundary_face.size(); i++){
        m1.push_back(boundary_face[i]);
        m2.push_back(boundary_face[i]);
        m1.push_back(1);
        m2.push_back(2);
        // cout<<boundary_face[i]<<" "<< nodes[l_face_to_node[m1]][0]<<" "<< nodes[l_face_to_node[m1]][1]<<endl;
        // cout<<boundary_face[i]<<" "<< nodes[l_face_to_node[m2]][0]<<" "<< nodes[l_face_to_node[m2]][1]<<endl;
        m1.clear();
        m2.clear();
    }
    for(int i=0; i<face_num; i++){
        P_faces.push_back(0);
    }
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
    }
    for(int i=0; i<center_num; i++){
        P_centers_old[i]=0;
        P_centers_new[i]=0;
    } 
    set_boundary_conditions();
    cout<<"5"<<endl;
    fvm();
    // cout<<l_face_to_cell.count();
    // v1_.push_back(9);
    // v1_.push_back(3);
    // cout<<l_cell_to_face[v1_]<<endl;
    // cout<<l_cell_to_face.size()<<endl;
    // for(int i=0; i<node_to_cell[50].size(); i++){
    //     cout<<node_to_cell[50][i]<<endl;
    // }
    // for(int i=1; i<node_num+1; i++){
    //     P_nodes[i-1] = cell_to_vertex_interpolation(i);
    // }
    ofstream MyFile("pressure_nodes.dat");
    MyFile << "VARIABLES= \"X\", \"Y\", \"K\"" <<endl;
    // for(int i=0; i<node_num; i++){
    //     double x = nodes[i+1][0]; 
    //     double y = nodes[i+1][1];
    //     // MyFile << x << " " << y <<" "<<(((-2)*sin(x+2*y))+(16*exp(2*x+3*y))-P_nodes[i])/3<<endl;
    //     MyFile << x << " " << y <<" "<<P_nodes[i]<<endl;
    //     // cout<<i+1<< " "<<P_faces[i]<<endl;
    // }
    double max = 0;
    for(int i=0; i<node_num; i++){
        MyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i]<<endl;
        if(abs(sin(nodes[i+1][0] + 2 * nodes[i+1][1]) + exp(2 * nodes[i+1][0] + 3 * nodes[i+1][1]) - P_nodes[i]) > max){
            max = abs(sin(nodes[i+1][0] + 2 * nodes[i+1][1]) + exp(2 * nodes[i+1][0] + 3 * nodes[i+1][1])  - P_nodes[i]);
        }
    }
    cout<<"max error: "<<max<<endl;
    // for(int i=0; i<face_num; i++){
    //     MyFile <<P_faces[i]<<endl;
    //     // cout<<i+1<< " "<<P_faces[i]<<endl;
    // }
    MyFile.close();
    ofstream MFile("pressure_analitycal.dat");
    MFile << "VARIABLES= \"X\", \"Y\", \"P\"" <<endl;
    for(int i=1; i<node_num+1; i++){
        MFile << nodes[i][0] << " " << nodes[i][1] <<" "<<sin(nodes[i][0] + 2 * nodes[i][1]) + exp(2 * nodes[i][0] + 3 * nodes[i][1])<<endl;
        // cout<<i+1<< " "<<P_faces[i]<<endl;
    }
    MFile.close();
    return 0;
}