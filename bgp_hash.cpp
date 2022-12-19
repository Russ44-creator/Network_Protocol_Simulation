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
#include <iterator>
#include <map>
#include <unordered_map>
#include <regex>
#include <sstream>
#include <set>
#include "method.h"

using namespace std;

#define INF INT_MAX //Infinity


// routing table
struct Routing{
	string name; //Router Name
	int N;  // Destination network N
	int d; // Distance
	string X;  //Next Router
    int as; 
};

const int sz=10001; //Maximum possible number of vertices. Preallocating space for DataStructures accordingly
vector<Routing> routers[100000];  // route table


vector<tuple<int,int,int>> ve[sz]; 
vector<tuple<int, int, int>> lsa_info[sz]; 

std::queue<int> ebgp_queue;
std::queue<int> ibgp_queue;

// enable the hash table use the pair as the key
template <class T> inline void hash_combine(size_t &seed, T const &v) {
    seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
struct pair_hash {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const {
        size_t seed = 0;
        hash_combine(seed, p.first);
        hash_combine(seed, p.second);
        return seed;
    }
};
std::unordered_map<pair<int, int>, int, pair_hash> bgp_map[sz];  // pair(this_port_num, next_port_num), distance
std::unordered_map<pair<int, int>, int, pair_hash> min_dist[sz]; // pair(next_port, target_edge_port), 3rd_int: cost
std::unordered_map<int, int> as_min_dist[sz];


void Dijkstra_ibgp(int source, int k);
void send_ibgp(int source, int end, int k);
bool send_ebgp(int source, int end, int k);


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
			rp.N = k / 2 * i + j;
			rp.d = 1;
			rp.X = '-';
            rp.as = -1;
			routers[i].push_back(rp);
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
    // edgeport
	for (int i = 0; i < layer_2; i++){
		for (int j = 0; j < pod/2; j++){
			// int distance = rand() % 10 + 1;
			int distance = 1;
			ve[i].push_back(make_tuple(layer_2 + j + (i / (k/2)) * (k/2), i, distance));
			ve[layer_2 + j + (i / (k/2)) * (k/2)].push_back(make_tuple(i, layer_2 + j + (i / (k/2)) * (k/2), distance));
		}		
	}
	// Convergence layer
	for (int i = layer_2; i < layer_2 * 2; i++){
		int group = (i - layer_2) % (k/2);
		for(int j = 0; j < k/2; j++){
			// int distance = rand() % 10 + 1;
			int distance = 1;
			ve[i].push_back(make_tuple(2*layer_2 + (k/2) * group + j, i, distance));
			ve[2 * layer_2 + (k/2) * group + j].push_back(make_tuple(i, 2 * layer_2 + (k/2) * group + j, distance));
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
	vector<Routing> r;
	for(int i = 0;i < k * k * 5 / 4;i++){
		r = routers[i];
		string temp = name_to_num(i, k);
		cout<<"***********************routing table "<<temp<<"************************\n";
		cout<<"To target N\t\tDistanced\t\tNext RouterX\t\tNext As\n";
		for(int i = 0;i < r.size();i++){
			cout<<r[i].N<<"\t\t\t"<<r[i].d<<"\t\t"<<r[i].as<<endl;
		}
	}
}


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


bool update_ebgp(int k){
	int edge = k*k/2;
	set<int> sent_port;
    bool hasEBGP = false;
    std::queue<int> queue_temp = ebgp_queue;
	while(!queue_temp.empty()){
		int source = ebgp_queue.front();
		ebgp_queue.pop();
        queue_temp.pop();
		set<int>::iterator iter;
		if((iter = sent_port.find(source)) != sent_port.end()){
			continue;
		}
		else{
			sent_port.insert(source);
		}
		if(source < k*k){
			int group =  (source - edge) % (k / 2);
			for(int i = 0; i < k / 2; i++){
				int core_target = edge * 2 + i + group * (k / 2);
				if (send_ebgp(source, core_target, k))
                    hasEBGP = true;
			}
		}
		else{
			int group = (source - k * k) / (k / 2);
			for(int i = 0; i < k; i++){
				int agg_target = edge + k / 2 * i + group;
				if (send_ebgp(source, agg_target, k))
                    hasEBGP = true;
			}
		}
	}
    return hasEBGP;
}


bool send_ebgp(int source, int end, int k){
	bool hasEdit = false;
	int edge = k * k / 2;
	int end_group = 0;
	if(source >= k*k){
		end_group = end - ((end-edge)/(k/2))*(k/2) - edge;
	}
	else{
		end_group = end - ((end-edge*2)/(k/2))*(k/2) - edge * 2;
	}
	std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter;
	iter = min_dist[source].begin();
	while(iter != min_dist[source].end()){
		std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter_end;
		iter_end = min_dist[end].begin();
		bool isExist = false;
		while(iter_end != min_dist[end].end()){
			if (iter_end->first.second == iter->first.second){
				int new_dist = iter->second + 1;
				
				if((iter->second + std::get<2>(ve[source][end_group])) < iter_end->second){
				//if((iter->second + 1) < iter_end->second){
					hasEdit = true;
					min_dist[end][make_pair(source, iter->first.second)] = iter->second + std::get<2>(ve[source][end_group]);
					min_dist[end].erase(iter_end);
				}
				isExist = true;
				break;
			}
			iter_end ++;
		}
		if(isExist == false){
			hasEdit = true;
			min_dist[end][make_pair(source, iter->first.second)] = iter->second + std::get<2>(ve[source][end_group]);
		}
		iter ++;
	}
	if (hasEdit == true){
		if(end >= k*k){
			ebgp_queue.push(end);
		}
		else if(end < k*k){
            ebgp_queue.push(end);
			// send ibgp first to agg layer
			int group = (end - edge) / (k / 2);
			for(int i = 0; i < k / 2; i++){
				int agg_target = edge + i + group * (k / 2);
				if(agg_target == end){
					continue;
				}
				else{
					send_ibgp(end, agg_target, k);
				}
			}
            // send ibgp to edge layer
            for(int i = 0;i < k / 2;i++){
                int edge_target = i + group * (k / 2);
                send_ibgp(end, edge_target, k);
            }
		}
	}
    return hasEdit;
}


//通常对IBGP邻居会加上next-hop-self
//通过IBGP学习到的路由，不能传递给其他的IBGP,是为了防止路由环路
//所有的IBGP router两两相连，组成一个full-mesh的网络。Full-mesh的连接数与节点的关系是n*(n-1)
void send_ibgp(int source, int end, int k){
	bool hasEdit = false;
	int edge = k * k / 2;
	std::unordered_map<pair<int,int>, int, pair_hash>::iterator iter;
	iter = min_dist[source].begin();
	while(iter != min_dist[source].end()){
		std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter_end;
		iter_end = min_dist[end].begin();
		bool isExist = false;
		while(iter_end != min_dist[end].end()){
			if (iter_end->first.second == iter->first.second){
				int new_dist = iter->second + as_min_dist[source][end];

				if((iter->second + as_min_dist[source][end]) < iter_end->second){
					hasEdit = true;
					min_dist[end][make_pair(source, iter->first.second)] = iter->second + as_min_dist[source][end];
					min_dist[end].erase(iter_end);
				}
				isExist = true;
				break;
			}		
			iter_end ++;
		}
		if(isExist == false){
			hasEdit = true;
			min_dist[end][make_pair(source, iter->first.second)] = iter->second + as_min_dist[source][end];
		}
		iter ++;
	}
	if (hasEdit == true && end >= k*k/2){
		ebgp_queue.push(end); // send ebgp
	}
}


void bgp(int k){
    int layer_2 = k * k / 2;
	int switch_in_pod = k;
	int core = (k / 2) ^ 2;
	int agg = k * k / 2;
	int edge = k * k / 2;
	int num_of_all = k * k * 5 / 4;
	// get the info of the neighbor from the initial info in each pod
	for(int i = 0; i < agg + edge; i++){
		 for(int j = 0; j < ve[i].size(); j++){
			if (get<0>(ve[i][j]) < k * k){
				bgp_map[i][make_pair(i, get<0>(ve[i][j]))] = std::get<2>(ve[i][j]);
			}
		 }
	}

	// loop k pods to get the min distance from the agg to the edge in the pod
	// IBGP
	for(int i = 0;i < edge; i++){
		for(int j = 0; j < k / 2; j++){
			bgp_map[edge + j + (i / (k/2)) * (k/2)] = send_lsa_hash(bgp_map[i], bgp_map[edge + j + (i / (k/2)) * (k/2)], i, edge + j + (i / (k/2))* (k/2));
		}
	}

	for(int i = 0;i < agg; i++){
		for(int j = 0; j < (k / 2); j++){
			bgp_map[i] = send_lsa_hash(bgp_map[edge + j + (i / (k/2))* (k/2)], bgp_map[i], edge + j + (i / (k/2))* (k/2), i);
		}
	}

	// write the ibgp to the lsdb
	// build_lsdb(k);

	// use dijkstra to calculate the shortest path in each pod
	for(int i = 0; i < agg + edge; i++){
		Dijkstra_ibgp(i, k);
	}
	
	// EBGP from agg to core
	// intimate getting the route from the initial 
	for(int i = edge; i < edge + agg; i++)
		for(int j = 0; j < ve[i].size(); j++){
			if (get<0>(ve[i][j]) >= k * k)
				bgp_map[i][make_pair(i, get<0>(ve[i][j]))] = std::get<2>(ve[i][j]);
		}

	for(int i = k^2; i < num_of_all; i++){
		for(int j = 0; j < ve[i].size(); j++)
			bgp_map[i][make_pair(i, get<0>(ve[i][j]))] = std::get<2>(ve[i][j]);
	}
	for(int i = k*k/2; i < k*k; i++){
		ebgp_queue.push(i);
	}
    while(true){
		if(update_ebgp(k)){ 
		}
        else{
            break;
        }
    }
	
    // IBGP to the edge layer
    for(int i = 0; i < edge; i++){
        int group = i / (k/2);
        for(int j = 0; j < k / 2; j++){
            int agg_port = j + (k / 2) * group;
            std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter;
	        iter = min_dist[agg_port].begin();
	        while(iter != min_dist[agg_port].end()){
                bool isExist = false;
                std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter_end;
		        iter_end = min_dist[i].begin();
                while(iter_end != min_dist[i].end()){
                    if (iter_end->first.second == iter->first.second){
                        if((iter->second + as_min_dist[agg_port][i]) < iter_end->second){
                            min_dist[i][make_pair(agg_port, iter->first.second)] = iter->second + as_min_dist[agg_port][i];
                            min_dist[i].erase(iter_end);
				        }
                        isExist = true;
                        break;
                    }
                    iter_end ++;
                }
                if(isExist == false){
                    min_dist[i][make_pair(agg_port, iter->first.second)] = iter->second + as_min_dist[agg_port][i];
		        }
                iter ++;
            }
        }
    }
}


int dis[sz]; //Stores shortest distance
int next_node[sz][sz];
bool vis[sz]={0}; //Determines whether the node has been visited or not


void Dijkstra_ibgp(int source, int k){
	bool vis[sz]={0};
	int server_num = k * k / 2;
	for(int i = 0;i < sz;i++) //Set initial distances to Infinity
        dis[i] = INF;
	 //Custom Comparator for Determining priority for priority queue (shortest edge comes first)
    class prioritize{public: bool operator ()(pair<int, int>&p1 ,pair<int, int>&p2){return p1.second>p2.second;}};
    priority_queue<pair<int,int> ,vector<pair<int,int>>, prioritize> pq; //Priority queue to store vertex,weight pairs
	pq.push(make_pair(source,dis[source]=0)); //Pushing the source with distance from itself as 0
	next_node[source][source] = source;
	while(!pq.empty()){
		pair<int, int> curr=pq.top(); //Current vertex. The shortest distance for this has been found
        pq.pop();
		int cv = curr.first,cw = curr.second; //'cw' the final shortest distance for this vertex
        if(vis[cv]) //If the vertex is already visited, no point in exploring adjacent vertices
            continue;
        vis[cv] = true;

		for(int i = 0;i < ve[cv].size();i++){ //Iterating through all adjacent vertices
			// destination: get<0>(ve[cv][i])
			if (get<0>(ve[cv][i]) < k*k){
				if(!vis[get<0>(ve[cv][i])] && get<2>(ve[cv][i]) + cw < dis[get<0>(ve[cv][i])]) {
				//If this node is not visited and the current parent node distance+distance from there to this node is shorter than the initial distace set to this node, update it
					pq.push(make_pair(get<0>(ve[cv][i]),(dis[get<0>(ve[cv][i])]=get<2>(ve[cv][i])+cw))); //Set the new distance and add to priority queue
					if (cv == source)
						next_node[source][get<0>(ve[cv][i])] = get<0>(ve[cv][i]);
					else
						next_node[source][get<0>(ve[cv][i])] = next_node[source][cv];
				}
			}
		}
	}
	int group = 0;
	if(source < k*k/2)
		group = source / (k/2);
	else
		group = (source - k*k/2) / (k/2);
	for(int i = 0; i < k/2; i++){
		if(dis[group*k/2+i] != INF){
			min_dist[source][make_pair(next_node[source][group*k/2+i], group*k/2+i)] = dis[group*k/2+i];
		}
	}
  
	// update the distance inside the AS
    for(int i = 0; i < k/2; i++)
        if(dis[k*k/2+group*k/2+i] != INF)
            as_min_dist[source][k*k/2+group*(k/2)+i] = dis[k*k/2+group*k/2+i];
}


void check_bgp(int k){
	int count = 0;
	for(int i = 0;i < k * k * 5 / 4;i++){		
		count += 1;
		if(count == k*k/2+1)
			cout << "edge++++++++++++" << endl;
		if(count == k*k+1)
			cout << "agg+++++++++++++" << endl;
		cout << min_dist[i].size()  << endl;
	}
}


int main(){
	int k = 20;
	// fattree(k);
	
	fattree_graph(k);
	
	int T=0;
	int count = 0;
	int switch_num = k * k * 5 / 4;
	int layer_2 = k * k / 2;
	auto start = std::chrono::system_clock::now();

	bgp(k);
	check_bgp(k);
    
	cout << "++++++++++++" << endl;
	int i = 0;
	std::unordered_map<pair<int,int> , int, pair_hash>::iterator iter;
	iter = min_dist[201].begin();
	while(iter != min_dist[201].end()){
		cout << iter->first.first<<" " <<iter->first.second<<" "<<iter->second << endl;
		i ++;
		iter++;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_spent = end-start;
	
	//printR2(k);
	cout << "\n time spent in OSPF: " << time_spent.count() << endl;;
	return 0;
}
