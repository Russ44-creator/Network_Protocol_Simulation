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

using namespace std;

#define INF INT_MAX //Infinity
const int sz=10001; // Maximum possible number of vertices. Preallocating space for DataStructures accordingly


// make the hash table can use the vector as the key
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


int dis[sz]; //Stores shortest distance
int next_node[sz][sz];
bool vis[sz]={0}; //Determines whether the node has been visited or not

void Dijkstra(int source, int k, vector<tuple<int, int, int>>lsdb[]){
	bool vis[sz]={0};
	int server_num = k * k / 2;
	for(int i = 0;i < sz;i++) //Set initial distances to Infinity
        dis[i] = INF;
	 //Custom Comparator for Determining priority for priority queue (shortest edge comes first)
    class prioritize{public: bool operator ()(pair<int, int>&p1 ,pair<int, int>&p2){return p1.second>p2.second;}};
    priority_queue<pair<int,int> ,vector<pair<int,int>>, prioritize> pq; //Priority queue to store vertex,weight pairs
	pq.push(make_pair(source,dis[source]=0)); //Pushing the source with distance from itself as 0
	int switch_num = k * k * 5 / 4;
	next_node[source][source] = -2;
	while(!pq.empty()){
		pair<int, int> curr=pq.top(); //Current vertex. The shortest distance for this has been found
        pq.pop();
		int cv = curr.first,cw = curr.second; //'cw' the final shortest distance for this vertex
        if(vis[cv]) //If the vertex is already visited, no point in exploring adjacent vertices
            continue;
        vis[cv] = true; 
        for(int i = 0;i < lsdb[cv].size();i++) //Iterating through all adjacent vertices
            if(!vis[get<0>(lsdb[cv][i])] && get<2>(lsdb[cv][i]) + cw < dis[get<0>(lsdb[cv][i])]) {
			//If this node is not visited and the current parent node distance+distance from there to this node is shorter than the initial distace set to this node, update it
                pq.push(make_pair(get<0>(lsdb[cv][i]),(dis[get<0>(lsdb[cv][i])]=get<2>(lsdb[cv][i])+cw))); //Set the new distance and add to priority queue
				if (cv == source)
					next_node[source][get<0>(lsdb[cv][i])] = get<0>(lsdb[cv][i]);
				else
					next_node[source][get<0>(lsdb[cv][i])] = next_node[source][cv];
			}
	}
}

