// Triangle mesh generator 
//NODES numeration -1 
#include <bits/stdc++.h>  
using namespace std;

const int n = 4;
const int n_side = n+1;
const int node_num=n_side*n_side;
const double step = 5/double(n);
map<int, vector<double> > nodes;
const int center_num=(n_side-1)*(n_side-1)*2;
map<int, vector<double> > centers;
const int face_num = ((n_side-1)*(n_side)*2) + (center_num/2);

map<vector<int> , int> l_cell_to_face;
map<vector<int> , int> l_face_to_cell;
map<vector<int> , int> l_face_to_node;
map<vector<int> , int> node_to_cell;
map<vector<int> , double> sn;


vector<double> P_nodes(node_num);
vector<double> P_faces(face_num);
vector<double> P_centers(center_num);
vector<double> P_centers_copy(center_num);

const double e=0.000001;


int main(){
    // cout<<face_num<<" :face num"<<endl;
    // cout<<center_num<<" :center num"<<endl;
    // cout<<node_num<<" :node num"<<endl;


    //Nodes coordinates generator
    vector<double> v; 
    int cnt = 1;
    for(double i=0; i<n_side; i++){
        for(double j=0; j<n_side; j++){
            v.push_back(j*step);
            v.push_back(i*step);
            nodes[cnt] = v;
            v.clear();
            cnt++;
        }
    }

    //Nodes Outputing 

    // map<int, vector<double> >::iterator it_node=nodes.begin(); 
    // while(it_node!=nodes.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it_node->first<<"-->"<<it_node->second[0]<<" "<<it_node->second[1]<<endl; 
    //     it_node++; 
    // } 


    //Cell to face generator
    vector<int> l; 
    l.push_back(n*n+1); 
    l.push_back(n*(2*n+1)+1); 
    l.push_back(1); 
    vector<int> cell_to_face_temp; 
    // cout << "============" << endl; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            l.at(1) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = l[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        l.at(0) += 1; 
        l.at(1) += 1; 
        l.at(2) += 1; 
    } 
    l.at(0) = 1; 
    l.at(1) = n*(n+1)+1; 
    l.at(2) = n*(2*n+1)+2; 
    for(int i = n*n+1; i < 2*n*n+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            l.at(2) += 1; 
        } 
        cell_to_face_temp.push_back(i); 
        for(int j = 1; j < 4; j++){ 
            cell_to_face_temp.push_back(j); 
            l_cell_to_face[cell_to_face_temp] = l[j-1]; 
            cell_to_face_temp.erase(cell_to_face_temp.end()-1); 
        } 
        cell_to_face_temp.clear(); 
        l.at(0) += 1; 
        l.at(1) += 1; 
        l.at(2) += 1; 
    } 
     
     
    // map<vector<int> , int>::iterator it=l_cell_to_face.begin(); 
    // while(it!=l_cell_to_face.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it->first[0]<<" "<<it->first[1]<<"-->"<<it->second<<endl; 
    //     it++; 
    // } 


    //Face to cell generator
    vector<int> v1; 
    for(int i=1; i<n_side+1; i++){
        v1.push_back(i);
        v1.push_back(1);
        l_face_to_cell[v1]=i;
        v1.erase(v1.end()-1);
        v1.push_back(2);
        l_face_to_cell[v1]=i;
        v1.clear();
    }
    int stopper = (n_side*2) -1;
    int start = n_side-1;
    int cell1 = n_side-1;
    int cell2 = 1;
    int temp;
    for(int i=n_side+1; i<face_num-n+1; i++){
        if(i == (start + stopper)){
            if(stopper == ((n_side*2) -1)){
                start = i+1;
                temp = cell1;
                cell1=cell2; 
                cell2 = temp+1;
                stopper = n_side-1;
                v1.push_back(i);
                v1.push_back(1);
                l_face_to_cell[v1]=cell1;
                v1.erase(v1.end()-1);
                v1.push_back(2);
                l_face_to_cell[v1]=cell1;
                v1.clear();
                cell1 = cell1+1;
                cell2= cell1-(n_side-1);
            }else if (stopper == n_side-1){
                start = i-1;
                stopper = ((n_side*2) -1);
                v1.push_back(i);
                v1.push_back(1);
                l_face_to_cell[v1]=cell2;
                v1.erase(v1.end()-1);
                v1.push_back(2);
                l_face_to_cell[v1]=cell2;
                v1.clear();
                cell2 = cell2; 
                cell1 = cell1 -1;  
            }
            continue;
        }

        if(stopper == (n_side*2) -1){
            temp = cell1;
            cell1=cell2; 
            cell2 = temp+1;
            v1.push_back(i);
            v1.push_back(1);
            l_face_to_cell[v1]=cell1;
            v1.erase(v1.end()-1);
            v1.push_back(2);
            l_face_to_cell[v1]=cell2;
            v1.clear();
        }else if(stopper == n_side-1){
            v1.push_back(i);
            v1.push_back(1);
            l_face_to_cell[v1]=cell1;
            v1.erase(v1.end()-1);
            v1.push_back(2);
            l_face_to_cell[v1]=cell2;
            v1.clear();
            cell1++;
            cell2++;
        }

    }
    for(int i=face_num-n+1; i<face_num+1; i++){
        v1.push_back(i);
        v1.push_back(1);
        l_face_to_cell[v1]= cell2;
        v1.erase(v1.end()-1);
        v1.push_back(2);
        l_face_to_cell[v1]= cell2;
        cell2++;
        v1.clear();
    }
    // map<vector<int> , int>::iterator it=l_face_to_cell.begin();
    // while(it!=l_face_to_cell.end()){
    //     cout<<it->first[0]<<" "<<it->first[1]<<"-->"<<it->second<<endl;
    //     it++;
    // }

    //Face to node generator
    vector<int> l1; 
    l1.push_back(n+2); 
    l1.push_back(2); 
    vector<int> face_to_node_temp; 
    for(int i = 1; i < n * n + 1; i++){ 
        if (i % n == 1 and i > 1){ 
            l1.at(0) += 1; 
            l1.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
    l1.at(0) = 2; 
    l1.at(1) = 1; 
    for(int i = n*n+1; i < n*(2*n+1)+1; i++){ 
        if (i % n == 1 and i != n*n+1){ 
            l1.at(0) += 1; 
            l1.at(1) += 1; 
        } 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
    l1.at(0) = 1; 
    l1.at(1) = n+2; 
    for(int i =  n*(2*n+1)+1; i < n*(3*n+2)+1; i++){ 
        face_to_node_temp.push_back(i); 
        for(int j = 1; j < 3; j++){ 
            face_to_node_temp.push_back(j); 
            l_face_to_node[face_to_node_temp] = l1[j-1]; 
            face_to_node_temp.erase(face_to_node_temp.end()-1); 
        } 
        face_to_node_temp.clear(); 
        l1.at(0) += 1; 
        l1.at(1) += 1; 
    } 
     
     
     
    // map<vector<int> , int>::iterator it1=l_face_to_node.begin(); 
    // while(it1!=l_face_to_node.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it1->first[0]<<" "<<it1->first[1]<<"-->"<<it1->second<<endl; 
    //     it1++; 
    // } 
    // cout << "============" << endl;

    //Node to cell generator 
    //cout << "Node to cell" << endl;
    vector<int> l2; 
    l2.push_back(2); 
    l2.push_back(n*n+1); 
    l2.push_back(n+1); 
    l2.push_back(n*(n+1)+1); 
    l2.push_back(n+2); 
    l2.push_back(n*n+2); 
    vector<int> node_to_cell_temp; 
    for(int i = n+3; i < n*(n+1); i++){ 
        if (i % (n+1) == 1 or i % (n+1) == 0){ 
            continue; 
        }else{ 
            node_to_cell_temp.push_back(i); 
            if (i % (n+1) == 2 and i != n+3){ 
                l2.at(0) += 1; 
                l2.at(1) += 1; 
                l2.at(2) += 1; 
                l2.at(3) += 1; 
                l2.at(4) += 1; 
                l2.at(5) += 1; 
            } 
            for(int j = 1; j < 7; j++){ 
                node_to_cell_temp.push_back(j); 
                node_to_cell[node_to_cell_temp] = l2[j-1]; 
                node_to_cell_temp.erase(node_to_cell_temp.end()-1); 
                 
            } 
            node_to_cell_temp.clear(); 
            l2.at(0) += 1; 
            l2.at(1) += 1; 
            l2.at(2) += 1; 
            l2.at(3) += 1; 
            l2.at(4) += 1; 
            l2.at(5) += 1; 
        } 
    }
    // map<vector<int> , int>::iterator it2=node_to_cell.begin(); 
    // while(it2!=node_to_cell.end()){ 
    //     // pair<vector<int> , int> res=*it; 
    //     // cout<<res.first.first<<" "<<res.first.second<<" = "<<res.second<<endl;// the same --> cout<<(*it)first<<"--> "<<(*it)second<<endl; 
    //     cout<<it2->first[0]<<" "<<it2->first[1]<<"-->"<<it2->second<<endl; 
    //     it2++; 
    // } 



    //Center Coordinate generator 
    vector<int> v_;
    vector<int> v1_;
    vector<vector<double> > points_;
    // cout << "-------------------" << endl;
    
    for(int i = 1; i < center_num + 1; i++){
        v_.push_back(i);
        v_.push_back(1);
        int face1 = l_cell_to_face[v_];
        v_.erase(v_.end()-1);
        v_.push_back(2);
        int face2 = l_cell_to_face[v_];
        v_.erase(v_.end()-1);
        v_.push_back(3);
        int face3 = l_cell_to_face[v_];
        v1_.push_back(face1);
        v1_.push_back(1);
        vector<double> node1 = nodes[l_face_to_node[v1_]];
        points_.push_back(node1);
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        vector<double> node2 = nodes[l_face_to_node[v1_]];
        points_.push_back(node2);
        v1_.clear();
        v1_.push_back(face2);
        v1_.push_back(1);
        vector<double> node3 = nodes[l_face_to_node[v1_]];

        if(count(points_.begin(), points_.end(), node3) == 0){
            points_.push_back(node3);
        }
        v1_.erase(v1_.end()-1);
        v1_.push_back(2);
        vector<double> node4 = nodes[l_face_to_node[v1_]];
        if(count(points_.begin(), points_.end(), node4) == 0){
            points_.push_back(node4);
        }
        v.push_back((points_[0][0] + points_[1][0] + points_[2][0]) / 3);
        v.push_back((points_[0][1] + points_[1][1] + points_[2][1]) / 3);
        centers[i]=v;
        points_.clear();
        v1_.clear();
        v_.clear();
        v.clear();
    }
    //Centers Outputing 
    // map<int, vector<double> >::iterator it_center=centers.begin(); 
    // while(it_center!=centers.end()){ 
    //     cout<<it_center->first<<"-->"<<it_center->second[0]<<" "<<it_center->second[1]<<endl; 
    //     it_center++; 
    // } 


    return 0;
}