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
#include <queue>
#include <climits>
#include <unordered_map>
#include <regex>
#include <sstream>
#include <iterator>
#include <map>
#include "method.h"

using namespace std;

#define INF INT_MAX //Infinity




// routing table
struct Routing{
	string name; //Router Name
	int N;  // Destination network N
	int d; // Distance
	string X;  //Next Router
};

vector<Routing> routers[1000000];  // route table
vector<tuple<int, int, int>> ve[sz]; // next_port_num, this_port_num, distance
std::unordered_map<pair<int, int>, int, pair_hash> lsa_map[sz]; // pair(this_port_num, next_port_num), distance
vector<tuple<int, int, int>> lsdb[sz]; // next_port_num, this_port_num, distance

std::queue<int> ebgp_queue;
std::queue<int> ibgp_queue;
std::unordered_map<pair<int, int>, int, pair_hash> bgp_map[sz];  // pair(this_port_num, next_port_num), distance
std::unordered_map<pair<int, int>, int, pair_hash> min_dist[sz]; // pair(next_port, target_edge_port), 3rd_int: cost
std::unordered_map<int, int> as_min_dist[sz];

// init fattree
void fattree(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;
	int layer_2 = k * k / 2;
	Routing rp;
	int N = core + pod * pod;
	for (int i = 0; i < layer_2; i++){
		auto temp = std::to_string(i);
		rp.name = 'A' + temp;
		for (int j = 0; j < k / 2; j++){
			rp.N = 2 * i + j;
			rp.d = 1;
			rp.X = '-';
			routers[i].push_back(rp);
		}	
	}
}
	
// init the fattree graph
void fattree_graph(int k){
	int core = (k / 2) ^ 2;
	int pod = k;
	int layer_2 = k * k / 2;
	int switch_in_pod = k;
	int server_each_pod = k * k /4;

	// access layer
	for (int i = 0; i < layer_2; i++){
		for (int j = 0; j < pod/2; j++){
			int distance = rand() % 10 + 1;
			ve[i].push_back(make_tuple(layer_2 + j + (i / (k/2)) * (k/2), i, distance));
			ve[layer_2 + j + (i / (k/2)) * (k/2)].push_back(make_tuple(i, layer_2 + j + (i / (k/2)) * (k/2), distance));
		}		
	}

	// Convergence layer
	for (int i = layer_2; i < layer_2 * 2; i++){
		int group = (i - layer_2) % (k/2);
		for(int j = 0; j < k/2; j++){
			int distance = rand() % 10 + 1;
			ve[i].push_back(make_tuple(2*layer_2 + (k/2) * group + j, i, distance));
			ve[2 * layer_2 + (k/2) * group + j].push_back(make_tuple(i, 2 * layer_2 + (k/2) * group + j, distance));
		}		
	}
}

// convert the port name to port number
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


// convert the port number to port name
string num_to_name(int source, int k){
		string name_router;
		int layer_2 = k * k / 2;
		if(source < layer_2){
			auto temp = std::to_string(source);
			name_router = 'A' + temp;
		}
		else if(source >= layer_2 && source < layer_2 * 2){
			auto temp = std::to_string(source - layer_2);
			name_router = 'B' + temp;
		}
		else{
			auto temp = std::to_string(source - layer_2 * 2);
			name_router = 'C' + temp;
		}
		return name_router;
}


// print the routing table
void printR(int k){
	vector<Routing> r;
	for(int i = 0;i < k * k * 5 / 4;i++){
		r = routers[i];
		string temp = name_to_num(i, k);
		cout<<"***********************routing table "<<temp<<"************************\n";
		cout<<"To target N\t\tDistanced\t\tNext RouterX\n";
		for(int i = 0;i < r.size();i++){
			cout<<r[i].N<<"\t\t\t"<<r[i].d<<"\t\t\t"<<r[i].X<<endl;
		}
	}
}


// print the graph
void printG(int k){
	cout<<"**********************Network Graph********************\n";
	for(int i = 0;i < k * k * 5 / 4;i++){
		cout << i << " ";
		for(int j = 0;j < ve[i].size();j++){
			if (get<2>(ve[i][j]) != INF){
				cout << "router " << get<0>(ve[i][j]) << " distance: " << get<2>(ve[i][j]) << ", ";
				if(j < ve[i].size() - 1) cout<<" ";
			}
		}
		cout << endl;
	}
}


// build and initialize the router's own lsa
void init_lsdb(int k){
	int num_of_all = k * k * 5 / 4;
	for(int i = 0; i < num_of_all; i++){
		 for(int j = 0; j < ve[i].size(); j++){
			lsa_map[i][make_pair(i, get<0>(ve[i][j]))] = std::get<2>(ve[i][j]);
		 }
	}
}


// send the lsa with hash map
std::unordered_map<pair<int,int> , int, pair_hash> send_lsa_hash(std::unordered_map<pair<int,int> , int, pair_hash>start, std::unordered_map<pair<int,int> , int, pair_hash> destination, int from, int to){
	std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter;
    iter = start.begin();
	std::unordered_map<pair<int,int> , int, pair_hash>::iterator it;
    while(iter != start.end()) {
		int token = 0;
		auto temp_pair = make_pair(iter->first.second, iter->first.first);	
    	if (destination.find(iter->first) != destination.end()) {
			token = 1;
			iter ++;
			continue;
    	}
		
		if(destination.find(temp_pair) != destination.end()){
			token = 1;
			iter ++;
			continue;
		}
		if (token == 0){
            destination.insert({iter->first, iter->second});
		}
		iter++;
	}
	return destination;
}        


// flood the lsa
void flood_lsa(int k){
	int core = (k / 2) ^ 2;
	int agg = k * k / 2;
	int edge = k * k / 2;
	int num_of_all = k * k * 5 / 4;

	// edge -> agg
	// edge just has the routing between edge and agg
	for(int i = 0;i < edge; i++){
		for(int j = 0; j < k / 2; j++){
			lsa_map[edge + j + (i / (k/2)) * (k/2)] = send_lsa_hash(lsa_map[i], lsa_map[edge + j + (i / (k/2)) * (k/2)], i, edge + j + (i / (k/2))* (k/2));
		}
	}
	
	// agg -> core
	// agg has the agg-core and agg-edge
	for(int i = 0;i < agg; i++){
		for(int j = 0; j < k / 2; j++){
			unordered_map<pair<int,int> , int, pair_hash> temp;
			temp = send_lsa_hash(lsa_map[i + edge], lsa_map[edge * 2 + j + (i % (k/2)) * (k / 2)], i + edge, edge * 2 + j + (i % (k/2)) * (k / 2));
			lsa_map[edge * 2 + j + (i % (k/2)) * (k / 2)] = temp;
		}
	}  

	
	// core -> agg
	// need de-duplicate
	for(int i = 0;i < agg; i++){
		for(int j = 0; j < k / 2; j++){
			unordered_map<pair<int,int> , int, pair_hash> temp;
			temp = send_lsa_hash(lsa_map[edge * 2 + j + (i % (k/2)) * (k / 2)], lsa_map[i + edge], edge * 2 + j + (i % (k/2)) * (k / 2), i + edge);
			lsa_map[i + edge] = temp;
		}
	}
	
	// other agg -> edge
	// no need to de-duplicate	
	for(int i = 0;i < agg; i++){
		for(int j = 0; j < (k / 2); j++){
			unordered_map<pair<int,int> , int, pair_hash> temp;
			temp = send_lsa_hash(lsa_map[edge + j + (i / (k/2))* (k/2)], lsa_map[i], edge + j + (i / (k/2))* (k/2), i);
			lsa_map[i] = temp;
		}
	}
	
	// edge -> agg
	// just one edge is enough
	// need to de-duplicate
	for(int i = 0; i < k; i++){
		for(int j = 0; j < k / 2; j++){
			lsa_map[edge + i * (k/2) + j] = send_lsa_hash(lsa_map[i * (k / 2)], lsa_map[edge + i * (k/2) + j], i * (k / 2), edge + i * (k/2) + j);
		}		
	}
	
	// agg -> core
	// just k / 2 agg port
	// need to de-duplicate	
	for(int i = 0; i < k / 2; i++){
		for(int j = 0; j < k / 2; j++){
			lsa_map[edge * 2 + j + (i % (k/2)) * (k / 2)] = send_lsa_hash(lsa_map[edge + i], lsa_map[edge * 2 + j + (i % (k/2)) * (k / 2)], edge + i, edge * 2 + j + (i % (k/2)) * (k / 2));
		}
	}	
}


// transfer the info in lsa_map to ve
void build_lsdb(int k){
    for(int i = 0;i < k * k * 5 / 4;i++){
		std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter;
		iter = lsa_map[i].begin();
		while(iter != lsa_map[i].end()) {
			auto temp_pair = make_pair(iter->first.second, iter->first.first);	
			lsdb[i].push_back(make_tuple(iter->first.second, i, iter->second));
			iter++;
		}
	}
}


//build the lsdb after the end of the flooding lsa
void build_routing_table(int k){
	int switch_num = k * k * 5 / 4;
	for(int source = 0; source < switch_num; source++){
		Dijkstra(source, k, lsdb);		
		//Printing final shortest distances from source
		for(int i = 0;i < k * k / 2; i++) {
			for (int j = 0; j < k / 2; j++){
				Routing rp;
				rp.name = num_to_name(source, k);
				rp.N = 2 * i + j;
				rp.d = dis[i] != INF? dis[i] + 1 :-1;
				if (rp.d == 1)
					rp.X = '-';
				else{
					rp.X = num_to_name(next_node[source][i], k);
				}
				routers[source].push_back(rp);
			}
		}	
	}
}


int main(){
	int k = 20;
	int method_type = 1; // 0: RIP, 1: OSPF, 2: BGP
	// fattree(k);
	fattree_graph(k);
	int T=0;
	int count = 0;
	int switch_num = k * k * 5 / 4;
	int layer_2 = k * k / 2;
	
	auto start = std::chrono::system_clock::now();

	init_lsdb(k);
	if (method_type == 1){
		flood_lsa(k);
		build_lsdb(k);
		build_routing_table(k);
	}
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_spent = end-start;
	
	//printR(k);
	cout << "\n time spent in OSPF: " << time_spent.count() << endl;;
	// printR(k);
	return 0;
}

