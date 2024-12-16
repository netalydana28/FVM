#include <bits/stdc++.h>  
using namespace std;
const int node_num=2601;
map<int, vector<double> > nodes;
const int center_num=2500;
map<int, vector<double> > centers;
int face_num;
map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_cell_to_node;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
set<int> l_node_to_nodes[node_num];
set<int> l_center_to_centers[center_num];
map<int, vector<int> > node_to_cell;
map<vector<int> , double> sn;
vector<int> boundary_face;
vector<int> nodes_boundary;
double l1_nodes[node_num];
double l1_centers[center_num];
double l2_nodes[node_num];
double l2_centers[center_num];
double u_star_nodes[node_num];
double v_star_nodes[node_num];
double u_n[node_num];
double v_n[node_num];
double u_n_1[node_num];
double v_n_1[node_num];
double P_nodes[node_num];
double P_nodes_copy[node_num];
double f[node_num];
vector<double> P_faces;
double P_centers[center_num];
double P_centers_copy[center_num];
double gauss_seidel_coeff[center_num][center_num];
double gauss_seidel_b[center_num];
double coeff[node_num][node_num];
double b[node_num];
double u = 1;
double v = 0;
double d = 0.1;
double delt = 0.000000001;
double ro = 1;
double left_boundary = 1; 
double right_boundary = 2; 
double up_boundary = 2; 
double bottom_boundary = 1; 
vector<double> split_double(string line, char delimeter){
    string my_str = "";
    vector<double> v;
    for(int i; i<line.size(); i++){
        if(((line[i] == delimeter) or (i==line.size()-1))){
            if(line[i+1] == delimeter){
                continue;
            }
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
        if((line[i] == delimeter) or (i==line.size()-1)){
            if(line[i+1] == delimeter){
                continue;
            }
            if(i==line.size()-1){
                my_str+=line[i];  
            }
            if(my_str[0]== 'q'){
                my_str = "";
                continue;
            }
            char* temp = new char[my_str.length()+1];
            strcpy(temp, my_str.c_str());
            int to_insert = atof(temp);
            v.push_back(to_insert);
            my_str = "";
        }else{
            my_str+=line[i];  
        }
    }
    return v;
}
void nodes_coor(){
    string myText;
    ifstream NodesFile("little_coor.txt");
    vector<double> v; 
    int cnt = 1;
    while (getline(NodesFile, myText)) {
        v = split_double(myText, ' ');
        // cnt = int(v[v.begin()]);
        v.erase(v.begin());
        v.erase(v.end()-1);
        nodes[cnt] = v;
        v.clear();
        cnt++;
    }
    NodesFile.close();
}
void cell_to_node(){
    string myText;
    ifstream NodesFile("little_links.txt");
    vector<int> v; 
    vector<int> v1; 
    while (getline(NodesFile, myText)) {
        v = split_int(myText, ' ');
        int cell = v[0];
        v.erase(v.begin());
        v.erase(v.begin());
        for(int i=1; i<5; i++){
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
void center_coordinates_generator(){//Checked without output
    vector<int> v;
    vector<double> v1;
    double x =0;
    double y =0;
    // double z =0;
    for(int i = 1; i < center_num + 1; i++){
        v.push_back(i);
        for(int j = 1; j<5; j++){
            v.push_back(j);
            register int node = l_cell_to_node[v];
            x += nodes[node][0];
            y += nodes[node][1];
            // z += nodes[node][2];
            v.erase(v.end()-1);
        }
        v1.push_back(x/double(4));
        v1.push_back(y/double(4));
        // v1.push_back(z/double(4));
        centers[i]=v1;
        v.clear();
        v1.clear();
        x = 0; 
        y = 0;
        // z = 0;
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
        for(int j=1; j<5; j++){
            key.push_back(j);
            register int trg_node = l_cell_to_node[key];
            trg_nodes.push_back(trg_node);
            key.erase(key.end()-1);
        }
        key.clear();
        register int m;
        for(int j=1; j<5; j++){
            if(j==4){
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
        int node1;
        int node2;
        while(s_it != m_it->first.end()){
            key.push_back(f_cnt); 
            key.push_back(s_cnt);
            l_face_to_node[key] = *s_it;
            if(s_cnt == 1){
                node1 = *s_it;
            }else if(s_cnt == 2){
                node2 = *s_it;
            }
            s_cnt++;
            s_it++;
            key.clear();
            
        }
        l_node_to_nodes[node1-1].insert(node2);
        l_node_to_nodes[node2-1].insert(node1);
        key.push_back(f_cnt);
        for(int i=1; i<key1.size()+1; i++){
            c_t_f[key1[i-1]-1].push_back(f_cnt);
            key.push_back(i);
            l_face_to_cell[key] = key1[i-1]; 
            key.erase(key.end()-1);
            if(key1.size() == 1){
                boundary_face.push_back(f_cnt);
                key.push_back(i+1);
                key.erase(key.end()-1);
                l_face_to_cell[key] = key1[i-1]; 
            }else{
                // cout<<key1[0]<<" "<<key1[1]<<endl;

                l_center_to_centers[key1[0]-1].insert(key1[1]);
                l_center_to_centers[key1[1]-1].insert(key1[0]);
            }
            // cout<<"Face "<<f_cnt<<" "<<key1[0]<< " "<<key1[1]<<endl;
           
            
        }
        f_cnt++;
        m_it++;
        key1.clear();
        key.clear();
    }
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
            key.erase(key.end()-1);
        }
        key.clear();
    }
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
void gauss_seidel(){
    for(int j=0; j<100; j++){
        for(int i=0; i<center_num; i++){
            P_centers_copy[i] = gauss_seidel_b[i];
            for(int j=0; j<center_num; j++){
                if (j == i){
                    continue;
                }
                P_centers_copy[i] = P_centers_copy[i] - (gauss_seidel_coeff[i][j] * P_centers[j]);
                P_centers[i] = (P_centers_copy[i]/gauss_seidel_coeff[i][i]);
            }
        }
    }
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
        double P = P_centers[cell_ind - 1]; 
        nom += P/d; 
        denom += 1/d; 
    } 
    return nom/denom; 
}
void fvm_convection(){
    for(int i=1; i<center_num+1; i++){
        vector<int> key;
        key.push_back(i);
        double ac=0;
        double bc=0;
        vector<double> center_coor = centers[i];
        for(int j=1; j<5; j++){
            key.push_back(j);
            register int face_ind = l_cell_to_face[key];
            vector<int> key1;
            key1.push_back(face_ind);
            double deltx;
            double delty;
            double nx; 
            double ny;
            double valid_denom;
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
                valid_denom=0;
                if(deltx != 0){
                    valid_denom += (nx/deltx);
                }
                if(delty != 0){
                    valid_denom += (ny/delty);
                }
                // cout<<" nx "<< nx<< " ny "<< ny<< " deltx "<< deltx<< " delty "<<delty<<" valid denom "<< valid_denom<<endl;
                bc += (area(face_ind))*P_faces[face_ind-1]*((d*valid_denom) - u*nx - v*ny);
                ac += (area(face_ind))*((d*valid_denom) - u*nx - v*ny);
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
                valid_denom=0;
                if(deltx != 0){
                    valid_denom += (nx/deltx);
                }
                if(delty != 0){
                    valid_denom += (ny/delty);
                }
                // cout<<" nx "<< nx<< " ny "<< ny<< " deltx "<< deltx<< " delty "<<delty<<" valid denom "<< valid_denom<<endl;
                gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind)*((d*valid_denom) - u*nx - v*ny));
                ac += (area(face_ind))*((d*valid_denom) - u*nx - v*ny);
                key1.clear();
            }
            key.erase(key.end()-1);
        }
        gauss_seidel_coeff[i-1][i-1] = ac;
        if(ac < 0){
            cout<<"YES"<<endl;
        }
        gauss_seidel_b[i-1] = bc;
        ac = 0;
        bc = 0; 
        key.clear();
    }
    gauss_seidel();
}
void set_boundary_conditions_u(){
    for(int i=0; i<face_num; i++){
        P_faces.push_back(0);
    }
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
        P_nodes_copy[i] = 0;
        u_n[i] = 0;
    }
    vector<int> v; 
    for(int i=0; i<boundary_face.size(); i++){
        int b_f = boundary_face[i];
        v.push_back(b_f); 
        v.push_back(1);
        int n1 = l_face_to_node[v];
        double x1 = nodes[l_face_to_node[v]][0];
        double y1 = nodes[l_face_to_node[v]][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        int n2 = l_face_to_node[v];
        double x2 = nodes[l_face_to_node[v]][0];
        double y2 = nodes[l_face_to_node[v]][1]; 
        if((fabs(x1-left_boundary)==0) and (fabs(x2-left_boundary)==0) and (y1>=1.4) and (y1<=1.6) and (y2>=1.4) and (y2<=1.6)){
        // if((fabs(x1-right_boundary)==0) and (fabs(x2-right_boundary)==0) and (y1>1.4) and (y1<1.6) and (y2>1.4) and (y2<1.6)){
        // if((fabs(y1-left_boundary)==0) and (fabs(y2-left_boundary)==0) and (x1>1.4) and (x1<1.6) and (x2>1.4) and (x2<1.6)){
        // if((fabs(y1-right_boundary)==0) and (fabs(y2-right_boundary)==0) and (y1>1.4) and (y1<1.6) and (y2>1.4) and (y2<1.6)){
            nodes_boundary.push_back(n1);
            nodes_boundary.push_back(n2);
            P_nodes[n1-1] = 0;
            P_nodes_copy[n1-1] = 0;
            P_nodes[n2-1] = 0;
            P_nodes_copy[n2-1] = 0;
            P_faces[b_f-1] = 0;
            u_n[n1-1] = 1;
            u_n[n2-1] = 1;
        }
        v.clear();
    }
    for(int i=0; i<center_num; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;
    } 
}
void set_boundary_conditions_v(){
    nodes_boundary.clear();
    for(int i=0; i<face_num; i++){
        P_faces[i] = 0;
    }
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
        P_nodes_copy[i] = 0;
        f[i] = 0;
        v_n[i] = 0;
    }
    vector<int> v; 
    for(int i=0; i<boundary_face.size(); i++){
        int b_f = boundary_face[i];
        v.push_back(b_f); 
        v.push_back(1);
        int n1 = l_face_to_node[v];
        double x1 = nodes[l_face_to_node[v]][0];
        double y1 = nodes[l_face_to_node[v]][1]; 
        v.erase(v.end()-1);
        v.push_back(2);
        int n2 = l_face_to_node[v];
        double x2 = nodes[l_face_to_node[v]][0];
        double y2 = nodes[l_face_to_node[v]][1]; 
        if((fabs(x1-left_boundary)==0) and (fabs(x2-left_boundary)==0) and (y1>=1.4) and (y1<=1.6) and (y2>=1.4) and (y2<=1.6)){
        // if((fabs(x1-right_boundary)==0) and (fabs(x2-right_boundary)==0) and (y1>1.4) and (y1<1.6) and (y2>1.4) and (y2<1.6)){
        // if((fabs(y1-left_boundary)==0) and (fabs(y2-left_boundary)==0) and (x1>1.4) and (x1<1.6) and (x2>1.4) and (x2<1.6)){
        // if((fabs(y1-right_boundary)==0) and (fabs(y2-right_boundary)==0) and (y1>1.4) and (y1<1.6) and (y2>1.4) and (y2<1.6)){
            nodes_boundary.push_back(n1);
            nodes_boundary.push_back(n2);
            P_nodes[n1-1] = 0;
            P_nodes_copy[n1-1] = 0;
            P_nodes[n2-1] = 0;
            P_nodes_copy[n2-1] = 0;
            P_faces[b_f-1] = 0;
            v[n1-1] = 0;
            v[n2-1] = 0;
        }
        v.clear();
    }
    for(int i=0; i<center_num; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;
    } 
}
void l1_find(){
    double du_dx;
    double du_dy;
    double d_2_u;//Poisson 2nd derivative
    for(int i=0; i<node_num; i++){
        double deltx = 0.02;  
        double cnt_x = 0; 
        double cnt_y = 0; 
        double cnt_2 = 0; 
        set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
        // cout<<"uij "<<i+1<<endl;
        // cout << "coor1 " << nodes[i+1][0] << " " << nodes[i+1][1] << endl;
        while(s_it != l_node_to_nodes[i].end()){ 
            cnt_2 += u_n[*s_it-1];
            if((nodes[*s_it][0] - nodes[i+1][0]) > 0){
                // if(nodes[i+1][0] != 2){
                    cnt_x += u_n[*s_it-1];
                    cnt_x -= (u_n[i]); 
                // }
            }
            if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                if(nodes[i+1][0] != 2){
                    cnt_y += u_n[*s_it-1]; 
                    cnt_y -= (u_n[i]);
                }
            }
            s_it++;
        } 
        // cout<<endl;

        // cnt_x -= (u_n[i]); 
        // cnt_y -= (u_n[i]);
        cnt_2 -= 4*(u_n[i]);

        du_dx = cnt_x/deltx; 
        du_dy = cnt_y/deltx; 
        d_2_u = cnt_2/(pow(deltx,2));

        l1_nodes[i] = -u * du_dx - v * du_dy + d_2_u;
    }
}
void l2_find(){
    double dv_dx;
    double dv_dy;
    double d_2_v;//Poisson 2nd derivative
    for(int i=0; i<node_num; i++){
        double deltx = 0.02;  
        double cnt_x = 0; 
        double cnt_y = 0; 
        double cnt_2 = 0; 
        set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
        // cout<<"uij "<<i+1<<endl;
        // cout << "coor1 " << nodes[i+1][0] << " " << nodes[i+1][1] << endl;
        while(s_it != l_node_to_nodes[i].end()){ 
            cnt_2 += v_n[*s_it-1];
            if((nodes[*s_it][0] - nodes[i+1][0]) > 0){
                // if(nodes[i+1][0] != 2){
                    cnt_x += v_n[*s_it-1];
                    cnt_x -= (v_n[i]); 
                // }
            }
            if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                if(nodes[i+1][0] != 2){
                    cnt_y += v_n[*s_it-1]; 
                    cnt_y -= (v_n[i]);
                }
            }
            s_it++;
        } 
        // cout<<endl;
        cnt_2 -= 4*(v_n[i]);

        dv_dx = cnt_x/deltx; 
        dv_dy = cnt_y/deltx; 
        d_2_v = cnt_2/(pow(deltx,2));

        l2_nodes[i] = -u * dv_dx - v * dv_dy + d_2_v;
    }
}
void u_star(){
    for(int i =0; i<node_num; i++){
        double x = nodes[i+1][0];
        double y = nodes[i+1][1];
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i+1) != nodes_boundary.end()){
            u_star_nodes[i] = 1;
        // }else if ((fabs(x-1)==0) and (y>4.5) and (y<6.5)){
        }else if ((fabs(x-left_boundary)==0) or (fabs(x-right_boundary)==0) or (fabs(y-bottom_boundary)==0) or (fabs(y-up_boundary)==0)){
            u_star_nodes[i] = 0;
        }else{
            u_star_nodes[i] = (delt * l1_nodes[i]) + u_n[i];
        }
    }
    // for(int i=0; i<node_num; i++){
    //     u_n[i] = u_star_nodes[i];
    // }
}
void v_star(){
    for(int i =0; i<node_num; i++){
        double x = nodes[i+1][0];
        double y = nodes[i+1][1];
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i+1) != nodes_boundary.end()){
            v_star_nodes[i] = 0;
        }else if ((fabs(x-left_boundary)==0) or (fabs(x-right_boundary)==0) or (fabs(y-bottom_boundary)==0) or (fabs(y-up_boundary)==0)){
            v_star_nodes[i] = 0;
        }else{
            v_star_nodes[i] = (delt * l2_nodes[i]) + v_n[i];
        }
    }
    // for(int i=0; i<node_num; i++){
    //     v_n[i] = v_star_nodes[i];
    // }
}
void P_find(){ 
    //Poisson equation 
    for(int i=0; i<node_num; i++){ 
        double deltx = 0.02;  
        double cnt_u = 0; 
        double cnt_v = 0; 
        set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
        // cout<<"uij "<<i+1<<endl;
        // cout << "coor1 " << nodes[i+1][0] << " " << nodes[i+1][1] << endl;
        while(s_it != l_node_to_nodes[i].end()){ 
            
            if((nodes[*s_it][0] - nodes[i+1][0]) > 0){
                // cout <<  nodes[*s_it][0] - nodes[i+1][0] << endl;
                cnt_u += u_star_nodes[*s_it-1];  
                // cout<<"ui+1j "<<*s_it<<" ";
                // cout << "coor2 " << nodes[*s_it][0] << " " << nodes[*s_it][1] << endl;
            }
            if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                // cout << nodes[*s_it][1] - nodes[i+1][1] << endl;
                cnt_v += v_star_nodes[*s_it-1];    
                // cout<<"uij+1 "<<*s_it<<" ";
                // cout << "coor2 " << nodes[*s_it][0] << " " << nodes[*s_it][1] << endl;
            }
            s_it++;
        } 
        // cout<<endl;

        cnt_u -= (u_star_nodes[i]); 
        cnt_v -= (v_star_nodes[i]); 

        // u_star_nodes[i] = cnt_u/deltx; 
        // v_star_nodes[i] = cnt_v/deltx; 
        // P_centers[i] = (ro/delt) * (cnt); 
        P_nodes[i] = (ro/delt) * (cnt_u + cnt_v)/deltx; 
    } 
    for(int i=0; i<node_num; i++){
        double x = nodes[i+1][0];
        double y = nodes[i+1][1];
        int i_1;
        if(fabs(x-left_boundary) == 0){ 
            P_nodes[i] = 0;
        }else if(fabs(x-right_boundary) == 0){
            // P_nodes[i] = 0;
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][0] - nodes[i+1][0]) < 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            // cout<<i<<" "<<i_1<<endl;
            // cout<<nodes[i+1][0]<<" "<<nodes[i_1+1][0]<<endl;
            P_nodes[i] = P_nodes[i_1];
        }else if(fabs(y-bottom_boundary) == 0){
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            // cout<<i<<" "<<i_1<<endl;
            // cout<<nodes[i+1][1]<<" "<<nodes[i_1+1][1]<<endl;
            P_nodes[i] = P_nodes[i_1];
        }else if(fabs(y-up_boundary) == 0){
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][1] - nodes[i+1][1]) < 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            P_nodes[i] = P_nodes[i_1];
            // cout<<i<<" "<<i_1<<endl;
            // cout<<nodes[i+1][0]<<" "<<nodes[i_1+1][0]<<endl;
        }
    }
    for(int i=0; i<node_num; i++){
        P_nodes_copy[i] = P_nodes[i];
    }
    ofstream UyFile("uv_derivative.dat");
    UyFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
    for(int i=0; i<node_num; i++){
        UyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes_copy[i] <<endl;
    }
    UyFile.close();
    // double to_check = 1;
    // for(int i=0; i<node_num; i++){
    //     P_nodes_copy[i] = P_nodes[i];
    // }
    double deltx = 0.02;
    for(int j=0; j<10000; j++){
        // to_check = find_max(P_nodes_copy); 
        // cout << "j " << j << endl;
        double cnt = 0;
        // int node_cnt = 0;
        for(int i=0; i<node_num; i++){
            if(l_node_to_nodes[i].size() != 4){
                // cout<<"node "<<i<<endl;
                // node_cnt++;
                continue;
            }
            set <int>::iterator s_it = l_node_to_nodes[i].begin();
            while(s_it != l_node_to_nodes[i].end()){
                cnt += P_nodes_copy[*(s_it)-1];
                s_it++;
            }
            P_nodes[i] = 0.25 * (cnt + f[i] * deltx * deltx);
            cnt = 0;
        }
        // cout<<"b_node "<<node_cnt<<endl;
        for(int i=0; i<node_num; i++){
            P_nodes_copy[i] = P_nodes[i];
        }
    }
    for(int i=0; i<node_num; i++){
        double x = nodes[i+1][0];
        double y = nodes[i+1][1];
        int i_1;
        if(fabs(x-left_boundary) == 0){ 
            P_nodes[i] = 0;
        }else if(fabs(x-right_boundary) == 0){
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][0] - nodes[i+1][0]) < 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            // cout<<nodes[i+1][0]<<" "<<nodes[i_1+1][0]<<endl;
            P_nodes[i] = P_nodes[i_1];
        }else if(fabs(y-bottom_boundary) == 0){
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            P_nodes[i] = P_nodes[i_1];
        }else if(fabs(y-up_boundary) == 0){
            set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
            while(s_it != l_node_to_nodes[i].end()){ 
                if((nodes[*s_it][1] - nodes[i+1][1]) < 0){
                    i_1 = *s_it-1; 
                }
                s_it++;
            }
            P_nodes[i] = P_nodes[i_1];
        }
    }
}
void uv_find(){
    double dP_dx;
    double dP_dy;
    for(int i=0; i<node_num; i++){
        double deltx = 0.02;  
        double cnt_u = 0; 
        double cnt_v = 0; 
        set <int>::iterator s_it = l_node_to_nodes[i].begin(); 
        while(s_it != l_node_to_nodes[i].end()){ 
            
            if((nodes[*s_it][0] - nodes[i+1][0]) > 0){
                cnt_u += P_nodes[*s_it-1];  
            }
            if((nodes[*s_it][1] - nodes[i+1][1]) > 0){
                cnt_v += P_nodes[*s_it-1];    
            }
            s_it++;
        } 
        // cout<<endl;

        cnt_u -= (P_nodes[i]); 
        cnt_v -= (P_nodes[i]); 

        dP_dx = cnt_u/deltx; 
        dP_dy = cnt_v/deltx; 

        u_n_1[i] = u_star_nodes[i] - ((delt/ro)*dP_dx);
        v_n_1[i] = v_star_nodes[i] - ((delt/ro)*dP_dy);
        // P_centers[i] = (ro/delt) * (cnt); 
        // P_nodes[i] = (ro/delt) * (u_star_nodes[i] + v_star_nodes[i])/deltx; 
    }
    for(int i=0; i<node_num; i++){
        double x = nodes[i+1][0];
        double y = nodes[i+1][1];
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i+1) != nodes_boundary.end()){ 
            u_n_1[i] = 1;
            v_n_1[i] = 0;
        }else if(fabs(x-left_boundary) == 0){
            u_n_1[i] = 0;
            v_n_1[i] = 0;
        }else if((fabs(x-right_boundary) == 0) or (fabs(y-bottom_boundary) == 0) or (fabs(y-up_boundary) == 0)){
            u_n_1[i] = u_star_nodes[i];
            v_n_1[i] = v_star_nodes[i];
        }
    }
    for(int i=0; i<node_num; i++){
        u_n[i] = u_n_1[i];
        v_n[i] = v_n_1[i];
    }
    ofstream UyFile("final_solution.dat");
    UyFile << "VARIABLES= \"X\", \"Y\", \"P1\", \"U1\", \"V1\"" <<endl;
    for(int i=0; i<node_num; i++){
        UyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i] << " " <<u_n_1[i]<<" "<<v_n_1[i]<<endl;
    }
}
int main(){
    nodes_coor();//collecting coordinates and links
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
    set_boundary_conditions_u();
    cout<<"5"<<endl;
    set_boundary_conditions_v();
    for(int m = 0; m<2; m++){
        l1_find();
        ofstream MFile("l1.dat");
        MFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
        for(int i=0; i<node_num; i++){
            MFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<l1_nodes[i]<<endl;
        }
        MFile.close();
        l2_find();
        ofstream MyFile("l2.dat");
        MyFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
        for(int i=0; i<node_num; i++){
            MyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<l2_nodes[i]<<endl;
        }
        MyFile.close();
        u_star();
        ofstream yFile("u_star.dat");
        yFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
        for(int i=0; i<node_num; i++){
            yFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<u_star_nodes[i]<<endl;
        }
        yFile.close();
        v_star();
        ofstream RyFile("v_star.dat");
        RyFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
        for(int i=0; i<node_num; i++){
            RyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<v_star_nodes[i]<<endl;
        }
        RyFile.close();
        P_find();
        ofstream PyFile("Poisson.dat");
        PyFile << "VARIABLES= \"X\", \"Y\", \"G\"" <<endl;
        for(int i=0; i<node_num; i++){
            PyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i]<<endl;
        }
        PyFile.close();
        uv_find();
    }

    // for(int i=1; i<node_num+1; i++){
    //     double x = nodes[i][0];
    //     double y = nodes[i][1];
    //     if((find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()) or (((fabs(x-left_boundary)==0) or (fabs(x-right_boundary)==0)) and (y>=bottom_boundary) and (y<=up_boundary))){
    //         // cout<<i<<endl;
    //         continue;
    //     }else{
    //         P_nodes[i-1] = cell_to_vertex_interpolation(i);
    //     }
    // }
    // for(int i=0; i<node_num; i++){
    //     l1_nodes[i] = P_nodes[i];
    // }
    // for(int i=0; i<center_num; i++){
    //     l1_centers[i] = P_centers[i];
    // }


    // fvm_convection();
    // for(int i=1; i<node_num+1; i++){
    //     double x = nodes[i][0];
    //     double y = nodes[i][1];
    //     if((find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()) or (((fabs(x-left_boundary)==0) or (fabs(x-right_boundary)==0)) and (y>=bottom_boundary) and (y<=up_boundary))){
    //         // cout<<i<<endl;
    //         continue;
    //     }else{
    //         P_nodes[i-1] = cell_to_vertex_interpolation(i);
    //     }
    // }
    // for(int i=0; i<node_num; i++){
    //     l2_nodes[i] = P_nodes[i];
    // }
    // for(int i=0; i<center_num; i++){
    //     l2_centers[i] = P_centers[i];
    // }


    // for(int i=1; i<node_num+1; i++){
    //     if(find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()){
    //         P_nodes[i-1] = 1;
    //     }else{
    //         P_nodes[i-1] = 0;
    //     }
    // }
    // u_star();
    // v_star();


    // P_find();

    // uv_find();
    return 0;
}