#include<iostream>
#include<fstream>
#include<algorithm>
#include<string>
#include<vector>
#include<cstdio>
#include<unistd.h>
#include<sys/time.h>
#include <ctime> 
#include <chrono>
#include <cstring>


using namespace std;
 
const int maxR = 6; // max number of router
 
 
// routing table
struct RIP{
	string name; //Router Name
	int N;  // Destination network N
	int d; // Distance
	string X;  //Next Router
};


// send package
vector<RIP> Send2(vector<RIP> r, string c){
	for(int i = 0;i < r.size();i++){
		r[i].d++;
		r[i].X = c;
	}
	return r;
}
 
vector<RIP> rip[100000];  // route table
char idToRname[] = "ABCDEF";
vector<int> ve[100000];  // network graph

// print the route table
void printR(){
	vector<RIP> r;
	for(int i = 0;i < maxR;i++){
		r=rip[i];
		cout<<"***********************routing table"<<idToRname[i]<<"************************\n";
		cout<<"To target N\t\tDistanced\t\tNext RouterX\n";
		for(int i=0;i<r.size();i++){
			cout<<r[i].N<<"\t\t\t"<<r[i].d<<"\t\t"<<r[i].X<<endl;
		}
	}
}

string name_to_num(int num, int k){
	int layer_2 = k * k / 2;
	if(num < layer_2){
			auto temp = std::to_string(num);
			return 'A' + temp;
		}
	else if(num >= layer_2 && num < layer_2 * 2){
		auto temp = std::to_string(num - layer_2);
		return 'B' + temp;
	}
	else{
		auto temp = std::to_string(num - layer_2 * 2);
		return 'C' + temp;
	}

}

void printR2(int k){
	vector<RIP> r;
	for(int i=0;i < k * k * 5 / 4;i++){
		r = rip[i];
		string temp = name_to_num(i, k);
		cout<<"***********************routing table "<<temp<<"************************\n";
		cout<<"To target N\t\tDistanced\t\tNext RouterX\n";
		for(int i = 0;i < r.size();i++){
			cout<<r[i].N<<"\t\t\t"<<r[i].d<<"\t\t"<<r[i].X<<endl;
		}
	}
}

void printG(int k){
	cout<<"**********************Network Graph********************\n";
	for(int i = 0;i < k * k * 5 / 4;i++){
		cout << i << " ";
		for(int j = 0;j < ve[i].size();j++){
			cout << ve[i][j];
			if(j < ve[i].size() - 1) cout<<" ";
		}
		cout << endl;
	}
}


// init fattree
void fattree(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;
	int layer_2 = k * k / 2;
	RIP rp;
	int N = core + pod * pod;
	for (int i = 0; i < layer_2; i++){
		auto temp = std::to_string(i);
		rp.name = 'A' + temp;
		for (int j = 0; j < k / 2; j++){
			rp.N = 2 * i + j;
			rp.d = 1;
			rp.X = '-';
			rip[i].push_back(rp);
		}
	}
}
	

void fattree_graph(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int layer_2 = k * k / 2;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;

	// access layer
	for (int i = 0; i < layer_2; i++){
		for (int j = 0; j < pod/2; j++){
			ve[i].push_back(layer_2 + j + (i / (k/2)) * (k/2));
		}		
	}
	// Convergence layer
	for (int i = layer_2; i < layer_2 * 2; i++){
		// access layer
		for(int j = 0; j < pod/2; j++){
			ve[i].push_back(j + (i - layer_2)/(k/2) * (k/2));
		}
		// core layer
		//for(int j = layer_2 * 2; j < layer_2*2 + core; j++){
		int group = (i - layer_2) % (k/2);
		for(int j = 0; j < k/2; j++){
			ve[i].push_back(2*layer_2 + (k/2) * group + j);
		}		
	}
	// core layer
	for(int i = layer_2 * 2; i < k * k * 5 / 4; i++){
		int group = (i - layer_2 * 2) % (k / 2);
		for(int j = 0; j < k; j++){
			ve[i].push_back(layer_2 + group + j * 2);
		}
	}
}

// init the routing table
void init(char* file){
	fstream infile(file);
	if(!infile.is_open()) cout<<"Fail to open file"<<endl;
	int N; // number of routers
	char name;
	RIP rp;
	infile >> N;
	for(int i = 0; i < N; i++){
		infile >> name;
		rp.name = name;
		while(infile>>rp.N>>rp.d>> rp.X&&rp.N ){  // 0: the end of one router
			rip[i].push_back(rp); // i: No. of the router
		}
	}
	infile.close();
}
 
// Initiate network graph
void initve(char* file){
	fstream infile(file);
	if(!infile.is_open()) cout<<"Fail to open file"<<endl;
	char v, v1;
	while(infile>>v) { 
		while(infile >> v1 && v1!='0'){
			ve[v-'A'].push_back(v1);
		}
	}
	infile.close();
}
 
bool cmp(RIP r1,RIP r2){
	return r1.N<r2.N;
}

bool update2(int k){
	bool hasEdit = false; 
	vector<RIP> send;
	vector<RIP> rip1[10000];
	int switch_num = k * k * 5 / 4;
	int layer_2 = k * k / 2;

	for(int i = 0;i < switch_num;i++) 
		rip1[i] = rip[i]; 

	for(int i = 0;i < switch_num;i++){
		string destination;
		if(i < layer_2){
			auto temp = std::to_string(i);
			destination = 'A' + temp;
		}
		else if(i >= layer_2 && i < layer_2 * 2){
			auto temp = std::to_string(i - layer_2);
			destination = 'B' + temp;
		}
		else{
			auto temp = std::to_string(i - layer_2 * 2);
			destination = 'C' + temp;
		}
		send = Send2(rip[i], destination); // send to neighbor
 
		for(int j = 0;j < ve[i].size();j++){	
			int lid = ve[i][j]; // neighbor No. 
			for(int x = 0;x < send.size();x++){
				int k;
				for(k = 0;k < rip1[lid].size();k++){
					if(send[x].N == rip1[lid][k].N){ // same destination
						if(send[x].X == rip1[lid][k].X){ // renew the distance if next jump is the same
							if(rip1[lid][k].d != send[x].d){
								rip1[lid][k] = send[x];
								hasEdit = true;
							} 
 
						}
						else{// shorter one
							if(send[x].d < rip1[lid][k].d){
								rip1[lid][k] = send[x];
								hasEdit = true;
							}
						}
						break;
					}
				}
 
				if(k == rip1[lid].size()){ // a new route
					rip1[lid].push_back(send[x]);
					hasEdit = true;
				}
			}			
		}
	}
    // write back
	for(int i = 0;i < switch_num;i++) {
		rip[i] = rip1[i];
		sort(rip[i].begin(),rip[i].end(),cmp);// sort with the destination
	}
	return hasEdit; 
}

 
int main(){
	int k = 4;
	fattree(k);
	//printR2(k);
	fattree_graph(k);
	//printG(k);
	int T=0;
	int count = 0;
	auto start = std::chrono::system_clock::now();
    
	while(true){
		sleep(10);
		cout<<"\n------------------------"<<++T<<" Renew------------------------\n";
		if(update2(k)){ 
			printR2(k);
		}
		else{
			auto end = std::chrono::system_clock::now();
			std::chrono::duration<double> time_spent = end-start;
			cout << "\n time spent in RIP: " << time_spent.count() << endl;;
			count += 1;
			cout<<"no renew\n";
			break;
			if (count >= 10){
				break; 
			}
		}
	}
	return 0;
}