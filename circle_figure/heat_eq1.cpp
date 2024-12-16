//Circle with triangle mesh
//Working with input data
#include <bits/stdc++.h>  
using namespace std;
const int node_num=5551;
map<int, vector<double> > nodes;
const int center_num=10520;
map<int, vector<double> > centers;
int face_num;
double delt = 0.2;
map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_cell_to_node;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<int, vector<int> > node_to_cell;
map<vector<int> , double> sn;
vector<int> boundary_face;
vector<int> nodes_boundary;
double k[center_num];
vector<double> k_boundary;
double P_nodes[node_num];
vector<double> P_faces;
double P_centers[center_num];
double P_centers_copy[center_num];
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
    ifstream NodesFile("obstacle_coor.txt");
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
    ifstream NodesFile("obstacle_links.txt");
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
        // cout<<node1<<" "<<node2<<" "<<node3<<endl;
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
                key.push_back(i+1);
                key.erase(key.end()-1);
                l_face_to_cell[key] = key1[i-1]; 
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
void gauss_seidel(){
    double to_check = 1;

    while(abs(to_check - find_max(P_centers))>e){
        // cout << find_max(P_centers) << endl;
        to_check = find_max(P_centers);
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
        // cout<<"HERE";
        // cout << abs(to_check - find_max(P_centers))<< endl;
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
void fvm(){
    for(int i=1; i<center_num+1; i++){
        vector<double> gDiff;
        vector<int> key;
        key.push_back(i);
        double ac=0;
        double bc=0;
        vector<double> center_coor = centers[i];
        for(int j=1; j<4; j++){
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
                bc += (area(face_ind))*P_faces[face_ind-1]*k_boundary[face_ind-1]*valid_denom;
                ac += (area(face_ind))*k_boundary[face_ind-1]*valid_denom;
                gDiff.push_back((-1)*(area(face_ind))*k_boundary[face_ind-1]*valid_denom);
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
                valid_denom=0;
                if(deltx != 0){
                    valid_denom += (nx/deltx);
                }
                if(delty != 0){
                    valid_denom += (ny/delty);
                }
                gauss_seidel_coeff[i-1][out_center_ind-1] = (-1)*(area(face_ind))*k[i-1]*valid_denom;
                ac += (area(face_ind))*k[i-1]*valid_denom;
                gDiff.push_back((-1)*(area(face_ind))*k[i-1]*valid_denom);

                key1.clear();

            }
            key.erase(key.end()-1);
        }
        gauss_seidel_coeff[i-1][i-1] = ac;
        gauss_seidel_b[i-1] = bc;
        ac = 0;
        bc = 0; 
        key.clear();
        gDiff.clear();
    }
    gauss_seidel();
}
void fvm2(){

}
void set_boundary_conditions(){
    for(int i=0; i<face_num; i++){
        P_faces.push_back(0);
    }
    for(int i=0; i<node_num; i++){
        P_nodes[i] = 0;
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
        if(((fabs(x1-1)<=0) and (fabs(x2-1)<=0))){
            nodes_boundary.push_back(n1);
            nodes_boundary.push_back(n2);
            P_nodes[n1-1] = 1;
            P_nodes[n2-1] = 1;
            P_faces[b_f-1] = 1;
        }
        v.clear();
    }
    for(int i=0; i<center_num; i++){
        P_centers[i]=0;
        P_centers_copy[i]=0;
    } 
    for(int i=0; i<center_num; i++){
        k[i] = 1;
    }
    for(int i=0; i<face_num; i++){
        k_boundary.push_back(1);
    }
}
int main(){
    nodes_coor();
    cout<<"1"<<endl;
    // cout<<nodes[3334][0]<<" "<<nodes[3334][1]<<endl;
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
    set_boundary_conditions();
    cout<<"5"<<endl;

    fvm();
    cout<<"6"<<endl;

    for(int i=1; i<node_num+1; i++){
        if(find(nodes_boundary.begin(), nodes_boundary.end(), i) != nodes_boundary.end()){
            // cout<<i<<endl;
            continue;
        }else{
            P_nodes[i-1] = cell_to_vertex_interpolation(i);
        }
    }
    ofstream MyFile("temperature_heat.dat");
    MyFile << "VARIABLES= \"X\", \"Y\", \"K\"" <<endl;
    for(int i=0; i<node_num; i++){
        MyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i]<<endl;
    }
    // for(int i=0; i<center_num; i++){
    //     MyFile << centers[i+1][0] << " " << centers[i+1][1] <<" "<<P_centers[i]<<endl;
    //     // cout<<i+1<< " "<<P_faces[i]<<endl;
    // }
    // for(int i=0; i<boundary_face.size(); i++){
    //     // MyFile <<P_faces[i]<<endl;
    //     cout<<i+1<< " "<<boundary_face[i]<<endl;
    // }
    // for(int i=0; i<node_num; i++){
    //     MyFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes[i]/delt<<endl;
    // }
    MyFile.close();
    // double P_nodes_old[node_num] = P_nodes;
    // fvm();
    // ofstream MFile("temperature_heat_time.dat");
    // MFile << "VARIABLES= \"X\", \"Y\", \"T\"" <<endl;
    // for(int i=0; i<node_num; i++){
    //     MFile << nodes[i+1][0] << " " << nodes[i+1][1] <<" "<<P_nodes_old[i]+(P_nodes[i]*delt)<<endl;
    // }
    // MFile.close();
    return 0;
}