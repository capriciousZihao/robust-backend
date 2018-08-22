// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira


#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <iostream>
#include <sstream>
#include <fstream> 
#include <string.h>
#include <vector>
#include <list>
#include <cstdlib>

#include <math.h>
#include <unistd.h>
#include <algorithm>  
#include <assert.h>

#include "g2o/core/eigen_types.h"
#include "g2o/types/slam2d/se2.h"
#include "g2o_interface_se2_zihao.hpp"

#include <Eigen/Dense>
#include "utils.hpp"

using namespace Eigen; 
using namespace std;
using namespace utils;

struct cluster
{
	int startLow, startHigh;
	int endLow, endHigh, nearestLCID=-1;
	int size, id_cluster;
	std::vector<std::array<double,6>> positionserial;
	std::vector<std::vector<std::pair<int, double> > > consistentGroup;
	std::vector<std::pair<int,int > > uncert_pair_group;
	std::vector< std::pair<double, std::pair<std::pair<int, int>, int> > > chiStatisWhenAdd2Cluster;
	// std::vector< std::pair<int,int>> nodesInCluster;

	cluster(): startLow(-1), startHigh(-1), endLow(-1), endHigh(-1), size(0) {}
	cluster(int start, int end, int idofcluster) : startLow(start), startHigh(start), endLow(end), endHigh(end), size(1), id_cluster(idofcluster){}

};

class Clusterizer
{
	typedef IntPairIDMap 			LoopToClusterIDMap;//loop is a pair of Vertex ID, this map is a int pair,LC, to a int, cluster ID 
	typedef IDintPairSetMap 		ClusterIDtoLoopsMap;
	double fiveNineBelief = 28, nineNineBelief = 44.8413, twelveBelief = 58.9413, upperLimitChi = 77.4, fiveupperLimitChi = 387, tenUpperLimit = 774,
		twoNineBelief = 11.345, threeNineBelief = 16.27, fourNineBelief = 21.11, PereentBelief95 = 7.8147, twoupperLimitChi = 154.8;
	double snrThres = 3.16; //3.16; 10
		// double snrThres = 10;

	std::vector<std::vector< std::vector<int> > >  node_sequence;
	std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > > transSequence_whole;
	std::vector<std::array<double,4>> VertexInf;

	std::map<std::pair<int, int>, std::array<double,9>> LC_Inf;	
	std::vector<double> distance;

	ClusterIDtoLoopsMap	clusterIDtoLoopsMap;
	
	std::vector<std::vector<int> > merge_vector;
	std::vector<int> killed_cluster;
	Matrix3d displayCov;
	int xx = 0;
	int wholeLoopN = 0;

public:
	int conflictNotInValid;
	std::vector<int>	vectorNewSplit,	vectorNewSplitCloser;
	LoopToClusterIDMap  loopToClusterIDMap, conflictClusterPairDistanceMap;
	std::vector<std::array<double,11>> OdoInf;
	std::vector<double> syntheticOdoAccumulate_dis;
	std::vector<std::pair<g2o::SE2, Matrix3d> > syntheticOdoAccumulate_trans_varMatrix;

	std::vector<std::pair<g2o::SE2, Matrix3d> > syntheticOdo10_trans_varMatrix, 
		syntheticOdo100_trans_varMatrix, syntheticOdo1000_trans_varMatrix, syntheticOdo10000_trans_varMatrix;
	std::vector<double > syntheticOdo10_dis, syntheticOdo100_dis, syntheticOdo1000_dis, syntheticOdo10000_dis;
	int NumOfRealLoops = -1;
	string     nameofclusterfile;
	std::vector<cluster> _clustersFound;
	std::set<std::pair<int,int> > set4GTL;
	std::vector<std::pair<int, std::vector<int> > > all_conflict_cluster;
	std::vector<std::pair<int, std::vector<int> > > all_uncertain_cluster;	
	std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	LP_Trans_Covar_Map;	
	find_element find_ele;
	// Assumes that in the vector the values are passed as (start_1,end_1), (start_2,end_2), ...
	int getClusterID(const IntPair& loop)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		LoopToClusterIDMap::iterator it;
		it = loopToClusterIDMap.find(loop);
		if(it == loopToClusterIDMap.end())
		{
			cout<<"the loop can not be found in map"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);	
		}
		// else
		// 	cout<<"check if get clUSter ID IS right "<<(*it).first.first<<" "<<(*it).first.second<<" id: "<<(*it).second <<endl;
		return loopToClusterIDMap[loop];
	}

	bool get_cluster_node_scequence(std::vector<std::vector< std::vector<int> > > & node_sequence, const std::vector<cluster> & clustersFound)
	{

		std::vector<int> nodes_dout, nodes_din;
		std::vector<std::vector<int> > element;
		bool retu = 0;
		node_sequence.clear();
		for(int i = 0; i < clustersFound.size(); ++i)
		{
			nodes_dout.clear();
			nodes_din.clear();	
			element.clear();	
			for(int j = 0; j < clustersFound[i].positionserial.size(); ++j)
			{
				nodes_dout.push_back(clustersFound[i].positionserial[j][0]);
				nodes_din.push_back(clustersFound[i].positionserial[j][3]);
			}
			element.push_back(nodes_dout);
			element.push_back(nodes_din);	
			node_sequence.push_back(element);
			// cout<<"i: "<<i<<endl;
			
			// cout<<"clustersFound[i].positionserial.size: "<<clustersFound[i].positionserial.size()<<endl;
			
			// cout<<"nodes_dout.size: "<<nodes_dout.size()<<endl;
			// cout<<"node_sequence[i][0].size: "<<node_sequence[i][0].size()<<endl;
			// cout<<"node_sequence[i][1].size: "<<node_sequence[i][1].size()<<endl;


			// std::cout << "nodes_dout sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_dout.begin(); it!=nodes_dout.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;

			// std::cout << "nodes_din  sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_din.begin(); it!=nodes_din.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;
		}
		// exit(0);
		retu = 1;
		return retu;
	}


	void setClusterID(const IntPair& loopClosure, int ID)
	{
		if(loopToClusterIDMap.find(loopClosure)!=loopToClusterIDMap.end())
		{
			int oldID = loopToClusterIDMap[loopClosure];

//			for(IntPairSet::iterator
//					it = clusterIDtoLoopsMap[oldID].begin(),
//					end = clusterIDtoLoopsMap[oldID].end();
//					it!=end;
//					it++)
//
//			{
				//if (*it == loopClosure){
					clusterIDtoLoopsMap[oldID].erase(loopClosure);
				//	loopToClusterIDMap[*it] = ID;
				//	break;
				//}
//			}
			clusterIDtoLoopsMap[ID].insert(loopClosure);
			loopToClusterIDMap[loopClosure] = ID;
			int clsuterSIZE = _clustersFound.size();
			if(ID > clsuterSIZE)
			{
				if(ID >_clustersFound.size() )
				{
					cout<<"cluster ID too big: "<<ID<<" while the whole clsuter number is "<<_clustersFound.size()<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
				}
			}
			// else
			// {
			// 	if(ID != -2)
			// 	{
			// 		cout<<"find the point that cluster ID is negative and not equal to -2: "<<ID<<endl;
			// 		printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			// 		exit(0);
			// 	}

			// }
		}
	}

	std::array<double,6> get_LC_Pos(int start, int end)
	{
		std::array<double,6> fullLoopInfo;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false;
		//get the position of start and end point std::vector<std::array<double,4>> VertexInf
		for(std::vector<std::array<double,4>>::const_iterator itVertex = VertexInf.begin(), 
			lendVertex = VertexInf.end();itVertex!=lendVertex;itVertex++)
		{
			if(int((*itVertex)[0])==start)
			{
                   fullLoopInfo[1]=(*itVertex)[1]; fullLoopInfo[2]=(*itVertex)[2]; //fullLoopInfo[2]=(*itVertex)[3];
				findStartVertexInMatrix=true;
			}
			if(int((*itVertex)[0])==end)
			{
				fullLoopInfo[4]=(*itVertex)[1]; fullLoopInfo[5]=(*itVertex)[2];// fullLoopInfo[2]=(*itVertex)[3];
				findEndVertexInMatrix=true;
			}
			if(findStartVertexInMatrix and findEndVertexInMatrix)
				break;
		}
		if(!(findStartVertexInMatrix and findEndVertexInMatrix))
			cout<<"can't find the position of the poses of the loop"<<endl;
		if(!findStartVertexInMatrix)
		{
			cout<<"start vertex id: "<<start<<endl;
		}
		if(!findEndVertexInMatrix)
		{
			cout<<"end vertex id: "<<end<<endl;
		}

		fullLoopInfo[0]=start;//  fullLoopInfo[1]=startPosition[0];  fullLoopInfo[2]=startPosition[1];
		fullLoopInfo[3]=end;  //  fullLoopInfo[4]=endPosition[0];    fullLoopInfo[5]=endPosition[1];
		return fullLoopInfo;
	}

	void prepare_to_signle_loop_pair_check_invalid_exist(std::pair<int, int> & loop_node_pair, IntPairSet::const_iterator & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2,
										   std::pair<int, int> & validVertexPoint)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;

		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (*LP_nodes).first;
		if(validVertexPoint.first > 0 and (abs(start - end) > validVertexPoint.first))
		{
			fultileBit1 = 1;
		}
		else
		{
			if(start > end)
			{
				// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
				accelerated_synthetic_odo( end, start, OdoInf, 
					result_synthetic_odo,  length[0], 
					syntheticOdo10_trans_varMatrix,
					syntheticOdo100_trans_varMatrix,
					syntheticOdo1000_trans_varMatrix,
					syntheticOdo10000_trans_varMatrix,
					syntheticOdo10_dis, syntheticOdo100_dis,
					syntheticOdo1000_dis, syntheticOdo10000_dis,
					fultileBit1);

				Trans1 = result_synthetic_odo.first.inverse();
				cov1   = result_synthetic_odo.second;
			}
			else
			{
				// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
				accelerated_synthetic_odo(start, end, OdoInf, 
					result_synthetic_odo,  length[0], 
					syntheticOdo10_trans_varMatrix,
					syntheticOdo100_trans_varMatrix,
					syntheticOdo1000_trans_varMatrix,
					syntheticOdo10000_trans_varMatrix,
					syntheticOdo10_dis, syntheticOdo100_dis,
					syntheticOdo1000_dis, syntheticOdo10000_dis,
					fultileBit1);

				Trans1 = result_synthetic_odo.first;
				cov1   = result_synthetic_odo.second;
			}	
		}
		//fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		if(fultileBit1 == 1)
		{
			if(validVertexPoint.first < 0)
				validVertexPoint.first = abs(start - end);
			else if (abs(start - end) < validVertexPoint.first)
			{
				validVertexPoint.first = abs(start - end);
			}
			return;
		}

		//synthesize the end odometry edges
		start = (*LP_nodes).second ;
		end   = loop_node_pair.second;

		if(validVertexPoint.second > 0 and (abs(start - end) > validVertexPoint.first))
		{
			fultileBit1 = 1;
		}
		else
		{
			if(start > end)
			{
				accelerated_synthetic_odo( end, start, OdoInf, 
					result_synthetic_odo,  length[2], 
					syntheticOdo10_trans_varMatrix,
					syntheticOdo100_trans_varMatrix,
					syntheticOdo1000_trans_varMatrix,
					syntheticOdo10000_trans_varMatrix,
					syntheticOdo10_dis, syntheticOdo100_dis,
					syntheticOdo1000_dis, syntheticOdo10000_dis,
					fultileBit2);	

				Trans2 = result_synthetic_odo.first.inverse();
				cov2   = result_synthetic_odo.second;
			}
			else
			{
				// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);
				accelerated_synthetic_odo(start, end, OdoInf, 
					result_synthetic_odo,  length[2], 
					syntheticOdo10_trans_varMatrix,
					syntheticOdo100_trans_varMatrix,
					syntheticOdo1000_trans_varMatrix,
					syntheticOdo10000_trans_varMatrix,
					syntheticOdo10_dis, syntheticOdo100_dis,
					syntheticOdo1000_dis, syntheticOdo10000_dis,
					fultileBit2);

				Trans2 = result_synthetic_odo.first;
				cov2   = result_synthetic_odo.second;
			}
		}
		if(fultileBit2 == 1)
		{
			if(validVertexPoint.second < 0)
				validVertexPoint.second = abs(start - end);
			else if (abs(start - end) < validVertexPoint.second)
			{
				validVertexPoint.second = abs(start - end);
			}
			return;
		}

		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(*LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans    = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];
	}

	void prepare_to_signle_loop_pair_check(std::pair<int, int> & loop_node_pair, IntPairSet::const_iterator & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;

		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (*LP_nodes).first;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[0], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit1);

			Trans1 = result_synthetic_odo.first.inverse();
			cov1   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[0], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit1);

			Trans1 = result_synthetic_odo.first;
			cov1   = result_synthetic_odo.second;
		}
		// //fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		// if(fultileBit == 1)
		// 	return;

		//synthesize the end odometry edges
		start = (*LP_nodes).second ;
		end   = loop_node_pair.second;

		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			// Trans2 = Trans2.inverse();

			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[2], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit2);

			// second_accelerated_synthetic_odo(end, start,  
			// 	result_synthetic_odo2, length[2], syntheticOdoAccumulate_trans_varMatrix,
			// 	syntheticOdoAccumulate_dis, fultileBit);		

			Trans2 = result_synthetic_odo.first.inverse();
			cov2   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);

			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[2], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit2);
			// second_accelerated_synthetic_odo(start, end,   
			// 	result_synthetic_odo2, length[2], syntheticOdoAccumulate_trans_varMatrix,
			// 	syntheticOdoAccumulate_dis, fultileBit);		

			Trans2 = result_synthetic_odo.first;
			cov2   = result_synthetic_odo.second;
		}

		// //fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		// if(fultileBit == 1)
		// 	return;

		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(*LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans    = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];

		// for(int i=0; i<4; i++){
		// 	cout<<" i: "<<i<<endl;
		// 	cout<<"cov: "<<endl;
		// 	cout<< FullInfo[i].second<<endl;
		// }
		// exit(0);
	}

	void prepare_to_signle_loop_pair_check(std::pair<int, int> & loop_node_pair, std::pair<int, int> & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit1,  bool & fultileBit2)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		bool fultileBit_;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1, Trans1_;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo, result_synthetic_odo2;
		Matrix3d  cov1, cov2, midCov, cov1_;
		bool fultileBit = 0;
		std::array<double ,5> length_;

		length[0][0] = 0; length[0][1] = 0; length[0][2] = 0;
		length[1][0] = 0; length[1][1] = 0; length[1][2] = 0;
		length[2][0] = 0; length[2][1] = 0; length[2][2] = 0;
		length[3][0] = 0; length[3][1] = 0; length[3][2] = 0;

		start = loop_node_pair.first;
		end   = (LP_nodes).first;
		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///
			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[0], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit1);

			Trans1 = result_synthetic_odo.first.inverse();
			cov1   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);
			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[0], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit1);

			Trans1 = result_synthetic_odo.first;
			cov1   = result_synthetic_odo.second;
		}
		//synthesize the end odometry edges
		start = (LP_nodes).second ;
		end   = loop_node_pair.second;

		if(start > end)
		{
			// synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			// Trans2 = Trans2.inverse();

			accelerated_synthetic_odo( end, start, OdoInf, 
				result_synthetic_odo,  length[2], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit2);

			Trans2 = result_synthetic_odo.first.inverse();
			cov2   = result_synthetic_odo.second;
		}
		else
		{
			// synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);

			accelerated_synthetic_odo(start, end, OdoInf, 
				result_synthetic_odo,  length[2], 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, syntheticOdo100_dis,
				syntheticOdo1000_dis, syntheticOdo10000_dis,
				fultileBit2);

			Trans2 = result_synthetic_odo.first;
			cov2   = result_synthetic_odo.second;
		}


		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans    = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);

		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];

	}

void prepare_to_signle_loop_pair_check_slow(std::pair<int, int> & loop_node_pair, IntPairSet::const_iterator & LP_nodes,
										   std::array<std::pair<g2o::SE2, Matrix3d>, 4> & FullInfo,
										   std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map,
										   std::array<std::array<double, 5>, 5> & length, bool & fultileBit)
	{
		//input variable: std::pair<int, int>  and full info of four segments.
		int start, end;
		g2o::SE2  Trans1, Trans2, midTrans, midTrans1;
		std::pair<g2o::SE2, Matrix3d > result_synthetic_odo;
		Matrix3d  cov1, cov2, midCov;
		fultileBit = 0;

		length[0][0] = 0;
		length[0][1] = 0;
		length[0][2] = 0;

		length[1][0] = 0;
		length[1][1] = 0;
		length[1][2] = 0;

		length[2][0] = 0;
		length[2][1] = 0;
		length[2][2] = 0;

		length[3][0] = 0;
		length[3][1] = 0;
		length[3][2] = 0;


		start = loop_node_pair.first;
		end   = (*LP_nodes).first;
		if(start > end)
		{
			synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1, length[0], fultileBit);///

			Trans1 = result_synthetic_odo.first.inverse();
		}
		else
		{
			synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1, length[0], fultileBit);

			Trans1 = result_synthetic_odo.first;
		}

		// //fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		// if(fultileBit == 1)
		// 	return;

		//synthesize the end odometry edges
		start = (*LP_nodes).second ;
		end   = loop_node_pair.second;

		if(start > end)
		{
			synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2, length[2], fultileBit);
			Trans2 = Trans2.inverse();
		}
		else
		{
			synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2, length[2], fultileBit);
		}

		// //fultileBit == 1 means the variance is larger then the distance, so we think there is no need to do next calculation
		// if(fultileBit == 1)
		// 	return;

		//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
		// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
		FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
		FullInfo[1] = LP_Trans_Covar_Map[(*LP_nodes)];
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
		midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
		midTrans = midTrans.inverse();
		midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
		FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

		length[1][0] = sqrt(FullInfo[1].first[0]*FullInfo[1].first[0]+FullInfo[1].first[1]*FullInfo[1].first[1]);
		length[1][1] = abs(FullInfo[1].first[0]);
		length[1][2] = abs(FullInfo[1].first[1]);
		length[1][3] = FullInfo[1].second(0,0);
		length[1][4] = FullInfo[1].second(1,1);


		length[3][0] = sqrt(FullInfo[3].first[0]*FullInfo[3].first[0]+FullInfo[3].first[1]*FullInfo[3].first[1]);
		length[3][1] = abs(FullInfo[3].first[0]);
		length[3][2] = abs(FullInfo[3].first[1]);
		length[3][3] = FullInfo[3].second(0,0);
		length[3][4] = FullInfo[3].second(1,1);
	

		length[4][0] = length[0][0]+length[1][0]+length[2][0]+length[3][0];
		length[4][1] = length[0][1]+length[1][1]+length[2][1]+length[3][1];
		length[4][2] = length[0][2]+length[1][2]+length[2][2]+length[3][2];

		// for(int i=0; i<4; i++){
		// 	cout<<" i: "<<i<<endl;
		// 	cout<<"cov: "<<endl;
		// 	cout<< FullInfo[i].second<<endl;

		// }
		// exit(0);

	}

	void find_cons_cluster(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;
			double bigChiErr = 0;
			if(_clustersFound[i].positionserial.size() > 1)
				limit = 2;

			//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
			double transformDistance = 100000000000.0;
			loop_node_pair.first = _clustersFound[i].positionserial.back()[0];
			loop_node_pair.second = _clustersFound[i].positionserial.back()[3];
			prepare_to_signle_loop_pair_check(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);


			if((futileBit1 or futileBit2) != 1)
			{	
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				bigChiErr = transformDistance;
				// if((reV.first == 1) and (futileBit == 0))
				if(reV.first == 1)
					pass  = 1;

				if(_clustersFound[i].positionserial.size() > 2)
				{
 					int si = _clustersFound[i].positionserial.size();
					loop_node_pair.first  = _clustersFound[i].positionserial[si-3][0];
					loop_node_pair.second = _clustersFound[i].positionserial[si-3][3];
					prepare_to_signle_loop_pair_check(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);
					reV_second_last = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual,
						length[4][0], futileBit);
					if(reV_second_last.second < transformDistance)
						transformDistance = reV_second_last.second;
					if(reV_second_last.second > bigChiErr)
						bigChiErr = reV_second_last.second;
					// if((reV_second_last.first == 1) and (futileBit == 0))
					if(reV_second_last.first == 1)
					{
						if(bigChiErr < upperLimitChi)
							pass  = 1;
					}
				}
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
					cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
					cout<<"transformDistance: "<<transformDistance<<endl;
					cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
					cout<<"cov:"<<endl;
					cout<<displayCov<<endl;
			}
			else
			{
				reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
				transformDistance = reV.second;
				cout<<" "<<endl;
				cout<<"fultileBit == 1"<<endl;
				cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
				cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
				cout<<"transformDistance: "<<transformDistance<<endl;
				cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				cout<<"cov:"<<endl;
				cout<<displayCov<<endl;

				if(transformDistance > fiveupperLimitChi )//tenUpperLimit
				{
					cout<<"add to conflict set although the NSR is higher than threshold"<<endl;
				}
			}
			//if current test cluster variance too big add it to the uncertain set
			if( abs(transformDistance - 100000000000.0) <0.001)
			{
				// uncertain_cluster.push_back(i);
				// continue;
			}
			//if pass chi2 test add this cluster id to the consistent set
			if(pass == 1)
			{
				cons_cluster_number.push_back(i);
				chiStatis.push_back(reV.second);			
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief twelveBelief
			else
			{
				if(transformDistance >  upperLimitChi )//tenUpperLimit
					conflict_cluster.push_back(i);
				// else if(transform_distance > 16.266)//%0.1
				// {
				// 	uncertain_cluster.insert(i);
				// }
				else 
				{
					// cons_cluster_number.push_back(i);
					// chiStatis.push_back(transform_distance);
					// reV.first = 1;
					// reV.second =transform_distance;
					// break;
					uncertain_cluster.push_back(i);
				}
			}

		}
	}
  
  
	void find_cons_cluster_second(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit,
		std::vector<std::pair<int,std::vector<int> > > & split, std::vector<std::pair<int,std::vector<double> > > & splitDis,
		std::map<int, double> & conflictClusterDistanceMap) // std::vector<std::vector<int> > & disForConsistent
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual, mind, maxd, maxdValid;
		bool futileBit1, futileBit2;
		std::vector<double> transDisWhole, distanceWholeValid;
		split.clear();
		splitDis.clear();
		std::pair<int,std::vector<int> > eleToSplit;
		std::pair<int, std::vector<double> > eleSplitDis;
		conflictClusterDistanceMap.clear();
		std::pair<int, int> validVertexPoint(-1,-1);

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;

			transDisWhole.clear();
			distanceWholeValid.clear();
			maxdValid = 0;

			for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[i].positionserial[wholeDis][0];
				loop_node_pair.second = _clustersFound[i].positionserial[wholeDis][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, 
					length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);

					if(transformDistance > maxdValid)
						maxdValid = transformDistance;
					// if((reV.first == 1) and (futileBit == 0))
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					distanceWholeValid.push_back(transformDistance);
					if(reV.first == 1)			
					{
						pass  = 1;
					}
				}
				else
				{
					cout<<" "<<endl;
					cout<<"fultileBit == 1"<<endl;
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);
					if(transformDistance > upperLimitChi )//tenUpperLimit
						cout<<"add to conflict set although the NSR is higher than threshold"<<endl;				
				}

			}

			mind = transDisWhole[0];//*(std::min_element(transDisWhole.begin(), transDisWhole.end()));
			// maxd = *(std::max_element(transDisWhole.begin(), transDisWhole.end()));
			maxd = 0;

			int minSerial = 0, maxSerial = 0;
			for(int findMax = 0; findMax < transDisWhole.size(); findMax++)
			{
				if(transDisWhole[findMax] > maxd)
				{
					maxd = transDisWhole[findMax];
					maxSerial = findMax;
				}
				if(transDisWhole[findMax] < transDisWhole[minSerial])
				{
					mind = transDisWhole[findMax];
					minSerial = findMax;
				}
			}

			cout<<"mind: "<<mind<<" "<<"maxd: "<<maxd<<endl;	
			cout<<"distanceWholeValid  serial: ";
			for(int toPrint = 0; toPrint < distanceWholeValid.size(); toPrint++)
			{
				cout<<distanceWholeValid[toPrint]<<" ";
				if(distanceWholeValid[toPrint] > maxdValid)
					maxdValid = distanceWholeValid[toPrint];
			}
			cout<<endl;
			cout<<"distance serial: ";
			for(int toPrint = 0; toPrint < transDisWhole.size(); toPrint++)
			{
				cout<<transDisWhole[toPrint]<<" ";
			}
			cout<<endl;

			//if pass chi2 test add this cluster id to the consistent set
			bool  futileBit69;
			std::vector<double> transformDistanceCluster, transformDistanceCluster2;
			if(pass == 1)
			{
				if(maxdValid > PereentBelief95)//fiveNineBelief
				{
					cout<<"clusterID: "<<i<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);
					cal_whole_cluster_distance_calculation(maxSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster, i);
					cout<<"distance to the biggest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						cout<<transformDistanceCluster[toPrint]<<" ";
					}
					cout<<endl;

					// transformDistanceCluster.clear();
					cal_whole_cluster_distance_calculation(minSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster2, i);
					cout<<"distance to the smallest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster2.size(); toPrint++)
					{
						cout<<transformDistanceCluster2[toPrint]<<" ";
					}
					cout<<endl;	
					eleToSplit.first = i;
					eleSplitDis.first = i;
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						if(transformDistanceCluster[toPrint] < transformDistanceCluster2[toPrint])
						{
							eleSplitDis.second.push_back(transformDistanceCluster[toPrint]);
							eleToSplit.second.push_back(toPrint);
						}
					}									
					split.push_back(eleToSplit);
					splitDis.push_back(eleSplitDis);
					// std::cin.get();
					// xx++;
				}
					cons_cluster_number.push_back(i);
					chiStatis.push_back(reV.second);		
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief
			else
			{
				// if(transformDistance >  fiveupperLimitChi )//tenUpperLimit  nineNineBelief   upperLimitChi twelveBelief
				if(maxd > upperLimitChi)
				{

					if(maxdValid < upperLimitChi and maxdValid > 0.00000000000000001)
					{
						conflictNotInValid++;
						cout<<"find conflict but the valid max is smaller than upper limit "<<endl;
						cout<<"current num in conflictNotInValid is "<<conflictNotInValid<<endl;
						cout<<"maxdValid: "<<maxdValid<<" maxd: "<<maxd<<endl;
						// conflictPairOnlyInWhole.push_back()

						std::cin.get();
					}
					conflict_cluster.push_back(i);
					if(maxd < upperLimitChi)
					{
						cout<<"the maxd is "<<maxd<<" which is samller than upperlimitchi"<<endl;
						// std::cin.get();
					}
					conflictClusterDistanceMap[i] = maxd;
				}
				else 
				{
					uncertain_cluster.push_back(i);
				}
			}

		}
		// wholeLoopN++;
	}

	void find_cons_cluster_second_invalid_exist(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit,
		std::vector<std::pair<int,std::vector<int> > > & split, std::vector<std::pair<int,std::vector<double> > > & splitDis,
		std::map<int, double> & conflictClusterDistanceMap) // std::vector<std::vector<int> > & disForConsistent
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual, mind, maxd, maxdValid;
		bool futileBit1, futileBit2;
		std::vector<double> transDisWhole, distanceWholeValid;
		split.clear();
		splitDis.clear();
		std::pair<int,std::vector<int> > eleToSplit;
		std::pair<int, std::vector<double> > eleSplitDis;

		conflictClusterDistanceMap.clear();

		std::pair<int, int> validVertexPoint(-1,-1);

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{

			bool pass  = 0;
			int limit = 1;

			transDisWhole.clear();
			distanceWholeValid.clear();
			maxdValid = 0;

			for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[i].positionserial[wholeDis][0];
				loop_node_pair.second = _clustersFound[i].positionserial[wholeDis][3];
				prepare_to_signle_loop_pair_check_invalid_exist(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, 
					length, futileBit1, futileBit2, validVertexPoint);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);

					if(transformDistance > maxdValid)
						maxdValid = transformDistance;
			
					// if((reV.first == 1) and (futileBit == 0))
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					distanceWholeValid.push_back(transformDistance);
					if(reV.first == 1)			
						pass  = 1;
				}
				else
					continue;				
			}

			mind = transDisWhole[0];//*(std::min_element(transDisWhole.begin(), transDisWhole.end()));
			// maxd = *(std::max_element(transDisWhole.begin(), transDisWhole.end()));
			maxd = 0;

			int minSerial = 0, maxSerial = 0;
			for(int findMax = 0; findMax < transDisWhole.size(); findMax++)
			{
				if(transDisWhole[findMax] > maxd)
				{
					maxd = transDisWhole[findMax];
					maxSerial = findMax;
				}
				if(transDisWhole[findMax] < transDisWhole[minSerial])
				{
					mind = transDisWhole[findMax];
					minSerial = findMax;
				}
			}

			// cout<<"distanceWholeValid  serial: ";
			for(int toPrint = 0; toPrint < distanceWholeValid.size(); toPrint++)
			{
				// cout<<distanceWholeValid[toPrint]<<" ";
				if(distanceWholeValid[toPrint] > maxdValid)
					maxdValid = distanceWholeValid[toPrint];
			}

			//if pass chi2 test add this cluster id to the consistent set
			bool  futileBit69;
			std::vector<double> transformDistanceCluster, transformDistanceCluster2;
			if(pass == 1)
			{
				if(maxdValid > PereentBelief95)//fiveNineBelief
				{
					cout<<"clusterID: "<<i<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);
					cal_whole_cluster_distance_calculation(maxSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster, i);

					// transformDistanceCluster.clear();
					cal_whole_cluster_distance_calculation(minSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster2, i);
						
					eleToSplit.first = i;
					eleSplitDis.first = i;
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						if(transformDistanceCluster[toPrint] < transformDistanceCluster2[toPrint])
						{
							eleSplitDis.second.push_back(transformDistanceCluster[toPrint]);
							eleToSplit.second.push_back(toPrint);
						}
					}									
					split.push_back(eleToSplit);
					splitDis.push_back(eleSplitDis);
				}
					cons_cluster_number.push_back(i);
					chiStatis.push_back(reV.second);		
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief
			else
			{
				// if(transformDistance >  fiveupperLimitChi )//tenUpperLimit  nineNineBelief   upperLimitChi twelveBelief

				// if(maxd > upperLimitChi)
				// {
				// 	if(maxdValid < upperLimitChi and maxdValid > 0.00000000000000001)
				// 	{
				// 		conflictNotInValid++;
				// 		cout<<"find conflict but the valid max is smaller than upper limit "<<endl;
				// 		cout<<"current num in conflictNotInValid is "<<conflictNotInValid<<endl;
				// 		cout<<"maxdValid: "<<maxdValid<<" maxd: "<<maxd<<endl;
				// 		// conflictPairOnlyInWhole.push_back()
				// 		std::cin.get();
				// 	}
				// 	conflict_cluster.push_back(i);
				// 	if(maxd < upperLimitChi)
				// 	{
				// 		cout<<"the maxd is "<<maxd<<" which is samller than upperlimitchi"<<endl;
				// 		// std::cin.get();
				// 	}
				// 	conflictClusterDistanceMap[i] = maxd;
				// }
				// else 
				// {
				// 	uncertain_cluster.push_back(i);
				// }
			}
		}
		// wholeLoopN++;
	}


	void find_cons_cluster_second_re(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::vector<int>  & conflict_cluster, std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, 
		std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug, std::vector<double> & chiStatis, bool & futileBit,
		std::vector<std::pair<int,std::vector<int> > > & split, std::vector<std::pair<int,std::vector<double> > > & splitDis,
		std::vector<int> & validMinDistance, std::map<int, double> & conflictClusterDistanceMap) // std::vector<std::vector<int> > & disForConsistent
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		chiStatis.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual, mind, maxd;
		bool futileBit1, futileBit2;
		std::vector<double> transDisWhole, distanceWholeValid;
		split.clear();
		splitDis.clear();
		std::pair<int,std::vector<int> > eleToSplit;
		std::pair<int, std::vector<double> > eleSplitDis;
		validMinDistance.clear();
		conflictClusterDistanceMap.clear();

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			bool pass  = 0;
			int limit = 1;

			transDisWhole.clear();
			distanceWholeValid.clear();

			for(int wholeDis = 0; wholeDis < _clustersFound[i].positionserial.size(); wholeDis++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[i].positionserial[wholeDis][0];
				loop_node_pair.second = _clustersFound[i].positionserial[wholeDis][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);

					distanceWholeValid.push_back(transformDistance);
					// if((reV.first == 1) and (futileBit == 0))
					cout<<" "<<endl;
					cout<<"fultileBit == 0"<<endl;
					if(reV.first == 1)			
					{
						pass  = 1;
						
					}
				}
				else
				{
					cout<<" "<<endl;
					cout<<"fultileBit == 1"<<endl;
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transDisWhole.push_back(transformDistance);
					if(transformDistance > upperLimitChi )//tenUpperLimit
						cout<<"add to conflict set although the NSR is higher than threshold"<<endl;				
				}

				// cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
				// cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
				// cout<<"transformDistance: "<<transformDistance<<endl;
				// cout<<"length: "<<length[0][0]<<" "<<length[1][0]<<" "<<length[2][0]<<" "<<length[3][0]<<" "<<length[4][0]<<" "<<endl;
				// 			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				// cout<<"cov:"<<endl;
				// cout<<displayCov<<endl;
			}

			mind = transDisWhole[0];//*(std::min_element(transDisWhole.begin(), transDisWhole.end()));
			// maxd = *(std::max_element(transDisWhole.begin(), transDisWhole.end()));
			maxd = 0;

			int minSerial = 0, maxSerial = 0;
			for(int findMax = 0; findMax < transDisWhole.size(); findMax++)
			{
				if(transDisWhole[findMax] > maxd)
				{
					maxd = transDisWhole[findMax];
					maxSerial = findMax;
				}
				if(transDisWhole[findMax] < transDisWhole[minSerial])
				{
					mind = transDisWhole[findMax];
					minSerial = findMax;
				}
			}
			double validMinDistanceELE = 100.0;
			cout<<"mind: "<<mind<<" "<<"maxd: "<<maxd<<endl;	
			cout<<"distanceWholeValid  serial: ";
			for(int toPrint = 0; toPrint < distanceWholeValid.size(); toPrint++)
			{
				cout<<distanceWholeValid[toPrint]<<" ";
				if(distanceWholeValid[toPrint] < validMinDistanceELE)
					validMinDistanceELE = distanceWholeValid[toPrint];
			}
			cout<<endl;
			cout<<"distance serial: ";
			for(int toPrint = 0; toPrint < transDisWhole.size(); toPrint++)
			{
				cout<<transDisWhole[toPrint]<<" ";
			}
			cout<<endl;

			//if pass chi2 test add this cluster id to the consistent set
			bool  futileBit69;
			std::vector<double> transformDistanceCluster, transformDistanceCluster2;
			if(pass == 1)
			{
				if(maxd > PereentBelief95)//fiveNineBelief
				{
					cout<<"clusterID: "<<i<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);
					cal_whole_cluster_distance_calculation(maxSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster, i);
					cout<<"distance to the biggest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						cout<<transformDistanceCluster[toPrint]<<" ";
					}
					cout<<endl;

					// transformDistanceCluster.clear();
					cal_whole_cluster_distance_calculation(minSerial,  _clustersFound,  
					LP_Trans_Covar_Map,  0, futileBit69, transformDistanceCluster2, i);
					cout<<"distance to the smallest error element: ";
					for(int toPrint = 0; toPrint < transformDistanceCluster2.size(); toPrint++)
					{
						cout<<transformDistanceCluster2[toPrint]<<" ";
					}
					cout<<endl;	
					eleToSplit.first = i;
					eleSplitDis.first = i;
					for(int toPrint = 0; toPrint < transformDistanceCluster.size(); toPrint++)
					{
						if(transformDistanceCluster[toPrint] < transformDistanceCluster2[toPrint])
						{
							eleSplitDis.second.push_back(transformDistanceCluster[toPrint]);
							eleToSplit.second.push_back(toPrint);
						}
					}									
					split.push_back(eleToSplit);
					splitDis.push_back(eleSplitDis);
					validMinDistance.push_back(validMinDistanceELE);
					// std::cin.get();
					// xx++;
				}
				else
				{
					cout<<"maxd: "<<maxd<<endl;
					cons_cluster_number.push_back(i);
					chiStatis.push_back(reV.second);
				}		
			}
			//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief
			else
			{
				// if(transformDistance >  fiveupperLimitChi )//tenUpperLimit nineNineBelief  upperLimitChi  twelveBelief
				if(maxd >  upperLimitChi )
				{
					conflict_cluster.push_back(i);
					if(maxd < upperLimitChi)
					{
						cout<<"the maxd is "<<maxd<<" which is samller than upperlimitchi"<<endl;
						// std::cin.get();
					}
					conflictClusterDistanceMap[i] = maxd;
				}
				else 
				{
					uncertain_cluster.push_back(i);
				}
			}

		}
		// wholeLoopN++;
	}


	void find_cons_cluster_whole_cluster_distance_calculation(
		IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, 
		std::vector<int> & cons_cluster_number,std::vector<int>  & conflict_cluster, 
		std::vector<int>  & uncertain_cluster, std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, 
		bool debug, std::vector<double> & chiStatis, bool & futileBit, 
		std::vector<double> & transformDistanceCluster, int & toTestClusterID)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();
		conflict_cluster.clear();
		uncertain_cluster.clear();
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;
		
		transformDistanceCluster.clear();

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		// for(int i=_clustersFound.size()-1; i >= 0; i--)
		int i = toTestClusterID;
		{
			bool pass  = 0;

			for(int wholeClusterTransformDistance = (_clustersFound[i].positionserial.size() -1); 
				wholeClusterTransformDistance >= 0; wholeClusterTransformDistance--)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[i].positionserial[wholeClusterTransformDistance][0];
				loop_node_pair.second = _clustersFound[i].positionserial[wholeClusterTransformDistance][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, LP_nodes, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(reV.first == 1 and (futileBit == 0))
					// if(reV.first == 1 )
						pass  = 1;
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 0"<<endl;
						cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}
				else
				{
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 1"<<endl;
						cout<<"test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}
				//if current test cluster variance too big add it to the uncertain set
				if( abs(transformDistance - 100000000000.0) <0.001)
				{
					uncertain_cluster.push_back(i);
					continue;
				}
				//if pass chi2 test add this cluster id to the consistent set
				if(pass == 1)
				{
					cons_cluster_number.push_back(i);
					chiStatis.push_back(reV.second);			
				}
				//if fail to pass the chi2 test, handle it based on the transform distance // nineNineBelief
				else
				{
					if(transformDistance >  tenUpperLimit)
						conflict_cluster.push_back(i);

					else 
					{
						uncertain_cluster.push_back(i);
					}
				}
			}
		}
	}

	void cal_whole_cluster_distance_calculation(
		int serial,  std::vector<cluster>& _clustersFound,  
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> > & LP_Trans_Covar_Map,  
		bool debug, bool & futileBit, 
		std::vector<double> & transformDistanceCluster, int & toTestClusterID)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		std::array<std::array<double, 5>, 5>  length;
		int start, end, eleSize = _clustersFound[toTestClusterID].positionserial.size();
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair, loop_node_pair_base; 
		std::pair<bool, double> reV, reV_second_last;
		double  transX_residual, transY_residual, transA_residual;
		bool futileBit1, futileBit2;
		
		transformDistanceCluster.clear();

		loop_node_pair_base.first = _clustersFound[toTestClusterID].positionserial[serial][0];
		loop_node_pair_base.second = _clustersFound[toTestClusterID].positionserial[serial][3];


		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		// for(int i=_clustersFound.size()-1; i >= 0; i--)

		{
			bool pass  = 0;
			for(int wholeClusterTransformDistance = 0; 
				wholeClusterTransformDistance < eleSize; wholeClusterTransformDistance++)
			{
				//calculate the transform distance, if there are more than two elements in the cluster, we will calculate the last two loop
				double transformDistance = 100000000000.0;
				loop_node_pair.first = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][0];
				loop_node_pair.second = _clustersFound[toTestClusterID].positionserial[wholeClusterTransformDistance][3];
				prepare_to_signle_loop_pair_check(loop_node_pair, loop_node_pair_base, FullInfo, LP_Trans_Covar_Map, length, futileBit1, futileBit2);

				if((futileBit1 or futileBit2) != 1)
				{	
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(reV.first == 1 and (futileBit == 0))
					// if(reV.first == 1 )
						pass  = 1;
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 0"<<endl;
						cout<<"test      loop: "<<(loop_node_pair_base).first<<" "<<(loop_node_pair_base).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}
				else
				{
					reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual, length[4][0], futileBit);
					transformDistance = reV.second;
					transformDistanceCluster.push_back(transformDistance);
					if(debug)
					{
						cout<<" "<<endl;
						cout<<"fultileBit == 1"<<endl;
						cout<<"test      loop: "<<(loop_node_pair_base).first<<" "<<(loop_node_pair_base).second<<endl;
						cout<<"loop : "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< toTestClusterID<<endl;
						cout<<"transformDistance: "<<transformDistance<<endl;
						cout<<"cov:"<<endl;
						cout<<displayCov<<endl;
					}
				}

			}
		}
	}

	bool judgeVarianceAndDistance(std::array<std::array<double, 5>, 5> & length)
	{
		bool accept = 1;
		for(int i =0; i < 4; i++)
		{
			for(int j =3; j<5; j++)
			{
				if(length[i][j] > length[i][0])
				{
					accept = 0;
					break;
				}	
			}
		}
		return accept;
	}

	bool judgeVarianceAndDistance_Small(std::array<double, 5> & length)
	{
		bool futileBitLocal = 0;

		// if(sqrt(sqrt(length[3]*length[3] + length[4]*length[4]) > 0.3*length[0]))
		// if(sqrt(length[3]*length[3] + length[4]*length[4]) > 13.25*length[0])

		// double SNR = (length[0]*length[0])/sqrt(length[3]*length[3] + length[4]*length[4]);
		double SNR = (length[0])/sqrt(length[3]*length[3] + length[4]*length[4]);

		// cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		if(SNR < snrThres)//3.16 corresponding to 5 dB
		// if(SNR < 1)
		// if(SNR < 10)
		{
			if(sqrt(length[3]*length[3] + length[4]*length[4]) > 0.1)
			{
				futileBitLocal = 1;
			}
		}

		return futileBitLocal;
	}

	void get_distance_to_neighbor_cluster(int i, std::vector<cluster>& _clustersFound, double & disx1, double & disy1, double & disx2, double & disy2)
	{
		disx1 = 10000;
		disy1 = 10000;
		disx2 = 10000;
		disy2 = 10000;
		for(int j = 0; j<_clustersFound.size(); j++)
		{
			if(j == i)
				continue;
			//disX1
			if(abs(_clustersFound[j].positionserial[0][1] - _clustersFound[i].positionserial.back()[1]) <  disx1)
				disx1 = abs(_clustersFound[j].positionserial[0][1] - _clustersFound[i].positionserial.back()[1]);
			if(abs(_clustersFound[j].positionserial.back()[1] - _clustersFound[i].positionserial.back()[1]) <  disx1)
				disx1 = abs(_clustersFound[j].positionserial.back()[1] - _clustersFound[i].positionserial.back()[1]);

			//disY1
			if(abs(_clustersFound[j].positionserial[0][2] - _clustersFound[i].positionserial.back()[2]) <  disy1)
				disy1 = abs(_clustersFound[j].positionserial[0][2] - _clustersFound[i].positionserial.back()[2]);
			if(abs(_clustersFound[j].positionserial.back()[2] - _clustersFound[i].positionserial.back()[2]) <  disy1)
				disy1 = abs(_clustersFound[j].positionserial.back()[2] - _clustersFound[i].positionserial.back()[2]);
			

			if(abs(_clustersFound[j].positionserial[0][4] - _clustersFound[i].positionserial.back()[4]) <  disx2)
				disx2 = abs(_clustersFound[j].positionserial[0][4] - _clustersFound[i].positionserial.back()[4]);
			if(abs(_clustersFound[j].positionserial.back()[4] - _clustersFound[i].positionserial.back()[4]) <  disx2)
				disx2 = abs(_clustersFound[j].positionserial.back()[4] - _clustersFound[i].positionserial.back()[4]);

			if(abs(_clustersFound[j].positionserial[0][5] - _clustersFound[i].positionserial.back()[5]) <  disy2)
				disy2 = abs(_clustersFound[j].positionserial[0][5] - _clustersFound[i].positionserial.back()[5]);	
			if(abs(_clustersFound[j].positionserial.back()[5] - _clustersFound[i].positionserial.back()[5]) <  disy2)
				disy2 = abs(_clustersFound[j].positionserial.back()[5] - _clustersFound[i].positionserial.back()[5]);
		}
	}
		
	std::array<double,2> chi2_test(const std::vector<double> & dis_clu, const std::vector<double> & dis_clu_real)
	{
		std::array<double,2> p_value = {0, 0};
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value[0] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		if(dis_clu.size()>2)
		{
			if(dis_clu_real.size() != dis_clu.size()-1)
			{
				cout<<"dis_clu_real.size: "<<dis_clu_real.size()<<endl;
				cout<<"dis_clu.size: "<<dis_clu.size()<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			
			sum = std::accumulate(std::begin(dis_clu_real), std::end(dis_clu_real), 0.0);  
			mean =  sum / dis_clu_real.size(); //均值  
			accum  = 0.0;  
			std::for_each (std::begin(dis_clu_real), std::end(dis_clu_real), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
			 chi_statis = (accum/mean); //卡方统计量
			p_value[1] = 1-chi2_p(dis_clu_real.size()-1, chi_statis);	
		}
		else
		{
			p_value[1] = p_value[0];
		}
		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		return p_value;
	}

	double chi2_test(const std::vector<double> & dis_clu)
	{
		double p_value ;
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value = 1-chi2_p(dis_clu.size()-1, chi_statis);	

		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		return p_value;
	}

	void updateSplit(std::vector<std::pair<int,std::vector<int> > > & split, 
		std::vector<std::pair<int,std::vector<double> > > & splitDis, std::vector<cluster> & _clustersFound,
		std::vector<int> & consCluster,
		std::vector<int> & conflictCluster,
		std::vector<std::pair<int,std::vector<int> > > & allConflictVector,
		std::vector<std::pair<int,std::vector<int> > > & allUncertainVector,
		int & consistentCluster)
	{
		if(split.size() > 1)
		{
			cout<<"split has "<<split.size()<< " consistent clusters, which is ";
			for(int toDisplay = 0; toDisplay < split.size(); toDisplay++)
			{
				cout<<split[toDisplay].first<<" ";
			}
			cout<<endl;
			cout<<"do split"<<endl;
			// std::cin.get();

			// if(split[toDisplay].first == )
			// {

			// }

			// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			// exit(0);
		}
		else
		{
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			return;
		}
		// assert(split.size() == 1);
		std::pair<double, std::pair<std::pair<int, int>, int> > ele_chiInfo;
		std::pair<int,int> changeLoop;
		std::pair<int,std::vector<int> > cofInfo;
		int conflictExist = find_ele.find(allConflictVector, consistentCluster), uncertainExist = find_ele.find(allUncertainVector, consistentCluster);


		if(consistentCluster != -1)
		{
			// cout<<"start to find abort"<<endl;
			// std::cin.get();
			for(int i = split[0].second.size() -1; i >= 0; i--)
			{
				//add position to cluster
				_clustersFound[consistentCluster].positionserial.push_back( _clustersFound[split[0].first].positionserial[split[0].second[i]]);
				//delete position serial
				if(consistentCluster == split[0].first)
				{
						printf("consistent clsuter equal to split cluster");
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
				}
				// _clustersFound[split[0].first].positionserial.erase(_clustersFound[split[0].first].positionserial.begin()+split[0].second[i]);

				//loop id change
				changeLoop.first  = _clustersFound[split[0].first].positionserial[split[0].second[i]][0];
				changeLoop.second = _clustersFound[split[0].first].positionserial[split[0].second[i]][3];
												if(consistentCluster > _clustersFound.size())
					{
						cout<<"find the point that cluster ID too big: "<<consistentCluster<<endl;
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}
				setClusterID(changeLoop, consistentCluster);

				//delete chi4
				// _clustersFound[split[0].first].chiStatisWhenAdd2Cluster.erase(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.begin()+split[0].second[i]);

				//add chi 
					ele_chiInfo.first = splitDis[0].second[i];
					if(abs(ele_chiInfo.first) < 0.000000000001)
						ele_chiInfo.first = std::max(splitDis[0].second[0], splitDis[0].second[1]);
					ele_chiInfo.second.first = changeLoop;
					ele_chiInfo.second.second = consistentCluster;

				_clustersFound[consistentCluster].chiStatisWhenAdd2Cluster.push_back(ele_chiInfo);				
			}

			for(int i = split[0].second.size() -1; i >= 0; i--)
			{
				//delete position serial
				if(i >1)
				{
					if(split[0].second[i] < split[0].second[i-1] )
					{
						printf("split[0].second[i] < split[0].second[i-1]");
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}
				}
				cout<<"split[0].second[i]: "<<split[0].second[i]<<endl;
				cout<<"positon size      : "<<_clustersFound[split[0].first].positionserial.size()<<endl;
				// std::cin.get();

				_clustersFound[split[0].first].positionserial.erase(_clustersFound[split[0].first].positionserial.begin()+split[0].second[i]);

				//delete chi4
				_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.erase(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.begin()+split[0].second[i]);	
			}
			// cout<<"end  to find abort"<<endl;
			// std::cin.get();

			//add conflict info 
			if(conflictCluster.size()>0)
			{
				if(conflictExist != -1)
				{
					if(allConflictVector[conflictExist].first != consistentCluster)
					{
						cout<<"consistentClusterID: "<<consistentCluster<<" allConflictVector[conflictExist].first"<<allConflictVector[conflictExist].first<<endl;
						printf("the two value should be equal.");
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}
					find_ele.add_vector_to_vector(conflictCluster, allConflictVector[conflictExist].second);
					if(uncertainExist != -1)
					{
						int haveConflict;
						for(int sizeCon = 0; sizeCon < conflictCluster.size(); sizeCon++)
						{
							// for(int toDelete = allUncertainVector[uncertainExist].second.size() - 1; toDelete >= 0; toDelete++)
							// {
							haveConflict = find_ele.find(allUncertainVector[uncertainExist].second, conflictCluster[sizeCon]);
							if(haveConflict != -1)
							{
								allUncertainVector[uncertainExist].second.erase(allUncertainVector[uncertainExist].second.begin()+haveConflict);
							}
							// }	
						}

					}
				}
				else
				{
					cofInfo.first = consistentCluster;
					cofInfo.second = conflictCluster;
					allConflictVector.push_back(cofInfo);
				}	
			}
			// if(consistentCluster ==95)
			// {
			// 	cout<<"consistentCluster is 97, so check if it introduce conflict cluster 95"<<endl;
			// 	// std::cin.get();
			// 	exit(0);
			// }
		}
		else
		{
			//new a cluster to hold the splited out loops
			cluster s(1,0,_clustersFound.size());
			
			for(int i = split[0].second.size() -1; i >= 0; i--)
			{
				s.positionserial.insert(s.positionserial.begin(), _clustersFound[split[0].first].positionserial[split[0].second[i]]);
				changeLoop.first  = _clustersFound[split[0].first].positionserial[split[0].second[i]][0];
				changeLoop.second = _clustersFound[split[0].first].positionserial[split[0].second[i]][3];
				//delete position serial
				_clustersFound[split[0].first].positionserial.erase(_clustersFound[split[0].first].positionserial.begin()+split[0].second[i]);
				//delete chi4
				_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.erase(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.begin()
					+split[0].second[i]);
				//delete loop from previous cluster	and delete the ID to the loop
				setClusterID(changeLoop, _clustersFound.size());
				// //add chi4 to new cluster
				// ele_chiInfo.first = splitDis[0].second[i];
				// ele_chiInfo.second.first = changeLoop;
				// ele_chiInfo.second.second = s.positionserial.size()-1;
				// s.chiStatisWhenAdd2Cluster.push_back(ele_chiInfo);
			}
			if(split[0].second.size() > 1)
			{
				for(int i = split[0].second.size() -1; i >= 0; i--)
				{
					//add chi4 to new cluster
					ele_chiInfo.first = splitDis[0].second[i];
					if(abs(ele_chiInfo.first) < 0.000000000001)
						ele_chiInfo.first = std::max(splitDis[0].second[0], splitDis[0].second[1]);
					ele_chiInfo.second.first = changeLoop;
					ele_chiInfo.second.second = i;
					s.chiStatisWhenAdd2Cluster.push_back(ele_chiInfo);
				}	
			}
		
			_clustersFound.push_back(s);
			//if there only one element remains, delete the chi4
			if(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.size() == 1)
			{
				_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.clear();
			}
			//add conflict info to the new split clsuter
			if(conflictCluster.size() > 0)
			{
				cofInfo.first  = _clustersFound.size()-1;
				cofInfo.second = conflictCluster;
				allConflictVector.push_back(cofInfo);
			}
		}

	}
void updateSplit_original(std::vector<std::pair<int,std::vector<int> > > & split, 
		std::vector<std::pair<int,std::vector<double> > > & splitDis, std::vector<cluster> & _clustersFound,
		std::vector<int> & consCluster,
		std::vector<int> & conflictCluster,
		std::vector<std::pair<int,std::vector<int> > > & allConflictVector,
		std::vector<std::pair<int,std::vector<int> > > & allUncertainVector,
		int & consistentCluster)
	{
		if(split.size() > 1)
		{
			cout<<"split has "<<split.size()<< " consistent clusters, which is ";
			for(int toDisplay = 0; toDisplay < split.size(); toDisplay++)
			{
				cout<<split[toDisplay].first<<" ";
			}
			cout<<endl;

			// if(split[toDisplay].first == )
			// {

			// }
			// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			// exit(0);
		}
		// assert(split.size() == 1);
		std::pair<double, std::pair<std::pair<int, int>, int> > ele_chiInfo;
		std::pair<int,int> changeLoop;
		std::pair<int,std::vector<int> > cofInfo;
		int conflictExist = find_ele.find(allConflictVector, consistentCluster);

		
			//new a cluster to hold the splited out loops
			cluster s(1,0,_clustersFound.size());
			
			for(int i = split[0].second.size() -1; i >= 0; i--)
			{
				s.positionserial.insert(s.positionserial.begin(), _clustersFound[split[0].first].positionserial[split[0].second[i]]);
				changeLoop.first  = _clustersFound[split[0].first].positionserial[split[0].second[i]][0];
				changeLoop.second = _clustersFound[split[0].first].positionserial[split[0].second[i]][3];
				//delete position serial
				_clustersFound[split[0].first].positionserial.erase(_clustersFound[split[0].first].positionserial.begin()+split[0].second[i]);
				//delete chi4
				_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.erase(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.begin()
					+split[0].second[i]);
				//delete loop from previous cluster	and delete the ID to the loop
				setClusterID(changeLoop, _clustersFound.size());
				// //add chi4 to new cluster
				// ele_chiInfo.first = splitDis[0].second[i];
				// ele_chiInfo.second.first = changeLoop;
				// ele_chiInfo.second.second = s.positionserial.size()-1;
				// s.chiStatisWhenAdd2Cluster.push_back(ele_chiInfo);
			}
			if(split[0].second.size() > 1)
			{
				for(int i = split[0].second.size() -1; i >= 0; i--)
				{
					//add chi4 to new cluster
					ele_chiInfo.first = splitDis[0].second[i];
					if(abs(ele_chiInfo.first) < 0.000000000001)
						ele_chiInfo.first = std::max(splitDis[0].second[0], splitDis[0].second[1]);
					ele_chiInfo.second.first = changeLoop;
					ele_chiInfo.second.second = i;
					s.chiStatisWhenAdd2Cluster.push_back(ele_chiInfo);
				}	
			}
		
			_clustersFound.push_back(s);
			//if there only one element remains, delete the chi4
			if(_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.size() == 1)
			{
				_clustersFound[split[0].first].chiStatisWhenAdd2Cluster.clear();
			}
			//add conflict info to the new split clsuter
			if(conflictCluster.size() > 0)
			{
				cofInfo.first  = _clustersFound.size()-1;
				cofInfo.second = conflictCluster;
				allConflictVector.push_back(cofInfo);
			}

	}
	void clusterize_zihao( const IntPairSet& loops, const char* filename)//, std::vector<int>& membership, int& clusterCount)
	{
		// std::map<int, std::set<int>> conflict, uncertain;
		std::map<int, double>  conflictClusterDistanceMap;
		double disSTART, disEND,  pre_p_value = 0.5;
		int num_loop=0, check = 0, consCluster_serialNum = -1;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false, ready2chiTest=false, futileBit = 0, fultileBit = 0;
		std::array<double,3> startPosition, endPosition, nearest_cluster;
		std::array<double,2> tem_dis, p_value;
		std::array<double,6> fullLoopInfo;//six elements:start ID,X,Y,end ID,X,Y
		std::pair<double, std::pair<std::pair<int, int>, int> > ele_chiInfo;
		std::set<int> splitSET, split2SET, spCONset, sp2CONset;

		std::vector<std::set<int> > conflict_cluster_set_vector;
		std::vector<int> conflict_cluster_set, uncertain_set,conflict_cluster_set_, uncertain_set_, 
			cons_cluster_number, cons_cluster_number2, cons_cluster_number_;
		std::pair<int,std::vector<int> > uncertain_conflict_ele;

		std::vector<double>  chiStatis, chiStatis2;
		std::vector<std::pair<int,std::vector<int> > >  split, split2 ;
		std::vector<std::pair<int,std::vector<double> > >  splitDis, splitDis2;

		//LC_Inf is a map, from nodes pair of loop closure to the six element array of 
		collect_vertexAndLP(filename, LC_Inf, LP_Trans_Covar_Map);

		prepare_accelerateBySynthetic(OdoInf, 
				syntheticOdo10_trans_varMatrix,
				syntheticOdo100_trans_varMatrix,
				syntheticOdo1000_trans_varMatrix,
				syntheticOdo10000_trans_varMatrix,
				syntheticOdo10_dis, 
				syntheticOdo100_dis,
				syntheticOdo1000_dis, 
				syntheticOdo10000_dis);

		// second_prepare_accelerated_synthetic_odo(OdoInf, syntheticOdoAccumulate_trans_varMatrix,
		// syntheticOdoAccumulate_dis);

		std::pair<g2o::SE2, Matrix3d>  result_synthetic_odo;

		std::array<double, 5>   length;

		if(loops.empty())
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			return;
		}
		_clustersFound.clear();

		// //debug check the result of one loop closure multiply the inverse of itself 
		// for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		// {

		// }
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{

			int start 	= it->first;
			int end 	= it->second;
			// if(start == 1575)
			// {
			// 	cout<<"find 1575: "<<endl;

			// 	exit(0);					
			// }

			// if(start == 2725)
			// 	exit(0);
			consCluster_serialNum = -2;
			if(start<end)
			{
				printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			//get loop closure vextexes position
			fullLoopInfo = get_LC_Pos( start,  end);
			//print loop number and cluster id of the loop
			num_loop++;
			cout<<"loop "<<num_loop<<" "<<start<<" "<<end<<endl;
			cout<<"current has "<<_clustersFound.size()<<" clusters"<<endl;

			if(_clustersFound.empty())
			{
				cluster s(start,end,_clustersFound.size());
				s.positionserial.push_back(fullLoopInfo);
				
				ele_chiInfo.first = 0;
				ele_chiInfo.second.first = *it;
				ele_chiInfo.second.second = s.positionserial.size()-1;
				s.chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//
				// s.nodesInCluster.push_back(*it);

				_clustersFound.push_back(s);

						cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
			" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;


				clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
				loopToClusterIDMap[*it] = _clustersFound.size()-1;


			}
			else
			{
				// if(it->first == 1441)
				// 	exit(0);
				//search for the nearest cluster to the loop
				cons_cluster_number.clear();
				chiStatis.clear();

				find_cons_cluster( it,  _clustersFound,  cons_cluster_number, conflict_cluster_set,
						uncertain_set, LP_Trans_Covar_Map, VertexInf, 1 , chiStatis,  futileBit);
				
				// find_cons_cluster_second( it,  _clustersFound,  cons_cluster_number, conflict_cluster_set,
				// 	uncertain_set, LP_Trans_Covar_Map, VertexInf, 1 , chiStatis,  futileBit, split, splitDis, conflictClusterDistanceMap);

				// find_cons_cluster_second_invalid_exist( it,  _clustersFound,  cons_cluster_number, conflict_cluster_set,
				// 	uncertain_set, LP_Trans_Covar_Map, VertexInf, 1 , chiStatis,  futileBit, split, splitDis, conflictClusterDistanceMap
				// 	);


				// assert(split.size() == 1);
				// assert(split.size() <= cons_cluster_number.size());
				cout<<"split size: "<<split.size()<<endl;
				cout<<"consistent size: "<<cons_cluster_number.size()<<endl;

				if(cons_cluster_number.size() == 0 and split.size() != 0)
				{
					cout<<"split size is not zero but cons cluster is zero: "<<endl;
					std::cin.get();
				}
				//display the conflict and uncertain set 
				cout<<"confilict cluster are: "<<endl;
				if(conflict_cluster_set.size() == 0)
					cout<<"empty"<<endl;
				else
				{
					for(int showConflict = 0; showConflict < conflict_cluster_set.size(); showConflict++)
					{
						cout<< conflict_cluster_set[showConflict]<<" ";
					}
				}
				cout<<endl;
				cout<<"uncertain cluster are: "<<endl;
				if(uncertain_set.size() == 0)
					cout<<"empty"<<endl;
				else
				{
					for(int showConflict = 0; showConflict < uncertain_set.size(); showConflict++)
					{
						cout<< uncertain_set[showConflict]<<" ";
					}
				}
				cout<<endl;
				// std::cin.get();
				int consClusterSerial = -1, uncertainClusterSerial = -1;
				bool exitBit = 0;
				std::vector<int > signToCalWholeCluster;

				//size equal to 0, then it means find no constent cluster, so construct a new cluster
				if(cons_cluster_number.size() == 0 )
				{
					cout<<"no consistent cluster found"<<endl;
					cluster s(start,end,_clustersFound.size());
					s.positionserial.push_back(fullLoopInfo);

					// s.chiStatisWhenAdd2Cluster.push_back( 0);
					// s.nodesInCluster.push_back(*it);
					ele_chiInfo.first = 0;
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = s.positionserial.size()-1;
					s.chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//

					_clustersFound.push_back(s);

					clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
					loopToClusterIDMap[*it] = _clustersFound.size()-1;
	

					consCluster_serialNum = _clustersFound.size()-1;
					// //construct the relationship between the clusters, whether it is confilict or uncertain
					find_ele.updateRelationship(consCluster_serialNum, conflict_cluster_set, 
						uncertain_set, cons_cluster_number, all_conflict_cluster, all_uncertain_cluster, exitBit, signToCalWholeCluster);

					// //add element to the cluster pair adn distance mapp
					// if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					// {
					// 	if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
					// 	{
					// 		printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// 		exit(0);		
					// 	}

					// 	for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
					// 	{
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
					// 	}	

					// }					
				}
				//find one constent cluster, add the loop to it
				else if (cons_cluster_number.size() == 1)
				{
					cout<<"only one consistent cluster: "<<cons_cluster_number[0]<<endl;

					consCluster_serialNum = cons_cluster_number[0];
					_clustersFound[consCluster_serialNum].positionserial.push_back(fullLoopInfo);

					clusterIDtoLoopsMap[consCluster_serialNum].insert(*it);
					loopToClusterIDMap[*it] = consCluster_serialNum;
					if(consCluster_serialNum > _clustersFound.size())
					{
						cout<<"find the point that cluster ID is too big: "<<consCluster_serialNum<<endl;
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}
		
					//save the chi statis when it pass the chi2 test and being add to this cluster, we save this for further check in intra check

					ele_chiInfo.first = chiStatis[0];
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = _clustersFound[consCluster_serialNum].positionserial.size()-1;
					_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//
					if(_clustersFound[consCluster_serialNum].positionserial.size() == 2)
						_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster[0].first = chiStatis[0];

					if(chiStatis.size() != 1)
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}

					find_ele.updateRelationship(consCluster_serialNum, conflict_cluster_set, 
						uncertain_set, cons_cluster_number, all_conflict_cluster, all_uncertain_cluster, exitBit, signToCalWholeCluster);

					// //add element to the cluster pair adn distance mapp
					// if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					// {
					// 	if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
					// 	{
					// 		printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// 		exit(0);		
					// 	}

					// 	for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
					// 	{
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
					// 	}	
					// }

				}
				//find more than one consistent cluster
				else 
				{
					cout<<"more than one consistent clusters"<<endl;
					int best_cons = 0;
					if(chiStatis.size() != cons_cluster_number.size())
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}

					cout<<" consistent cluster: "<<cons_cluster_number[0]<<endl;

					for(int selectFromMultiConsCluster = 1; selectFromMultiConsCluster < chiStatis.size(); selectFromMultiConsCluster++)
					{
						cout<<" consistent cluster: "<<cons_cluster_number[selectFromMultiConsCluster]<<endl;

						if(chiStatis[selectFromMultiConsCluster] < chiStatis[best_cons])
							best_cons = selectFromMultiConsCluster;
					}

					consCluster_serialNum = cons_cluster_number[best_cons];
					_clustersFound[consCluster_serialNum].positionserial.push_back(fullLoopInfo);
					cout<<"the final cluster selected is "<<consCluster_serialNum<<endl;

					clusterIDtoLoopsMap[consCluster_serialNum].insert(*it);
					loopToClusterIDMap[*it] = consCluster_serialNum;

					if(consCluster_serialNum > _clustersFound.size())
					{
						cout<<"find the point that cluster ID is too big: "<<consCluster_serialNum<<endl;
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}

					//save the chi statis when it pass the chi2 test and being add to this cluster, we save this for further check in intra check

					ele_chiInfo.first = chiStatis[best_cons];
					ele_chiInfo.second.first = *it;
					ele_chiInfo.second.second = _clustersFound[consCluster_serialNum].positionserial.size()-1;
					_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster.push_back( ele_chiInfo);//

					if(_clustersFound[consCluster_serialNum].positionserial.size() == 2)
						_clustersFound[consCluster_serialNum].chiStatisWhenAdd2Cluster[0].first = chiStatis[best_cons];

					//put other consistent cluster into uncertain set
					for(int selectFromMultiConsCluster = 0; selectFromMultiConsCluster < chiStatis.size(); selectFromMultiConsCluster++)
					{
						if(selectFromMultiConsCluster == best_cons)
							continue;
						uncertain_set.push_back(cons_cluster_number[selectFromMultiConsCluster]);
					}

					//update relationship

					find_ele.updateRelationship(consCluster_serialNum, conflict_cluster_set, 
						uncertain_set, cons_cluster_number, all_conflict_cluster, all_uncertain_cluster, exitBit, signToCalWholeCluster);

					// //add element to the cluster pair adn distance mapp
					// if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					// {
					// 	if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
					// 	{
					// 		printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					// 		exit(0);		
					// 	}

					// 	for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
					// 	{
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], consCluster_serialNum)] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
					// 		conflictClusterPairDistanceMap[std::pair<int,int> (consCluster_serialNum, conflict_cluster_set[addMap])] =
					// 		 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
					// 	}
					// }


					for(int piano = 0; piano < signToCalWholeCluster.size(); piano++)
					{
						int toTestClusterID = signToCalWholeCluster[piano];
						if(_clustersFound[toTestClusterID].positionserial.size() == 1)
						{
							if((_clustersFound.size() >  toTestClusterID) and (find_ele.find(cons_cluster_number, toTestClusterID) != -1))
							{
								int first_local  = _clustersFound[toTestClusterID].positionserial[0][0];
								int second_local = _clustersFound[toTestClusterID].positionserial[0][3];		
								IntPairSet::const_iterator it_ = loops.begin(), lend_ = loops.end();

								for(;it_!=lend_;it_++)
								{
									if(((*it_).first == first_local and (*it_).second == second_local) or 
										((*it_).first == second_local and (*it_).second == first_local))
										break;
								}
								if(it_==lend_)
								{
						    		printf("can not find the loop in cluster 29. This error is in %s on line %d\n",  __FILE__, __LINE__);
						    		exit(0);
								}

								cout<<"check cluster 29 each elemnt transform"<<endl;
								std::vector<double> transformDistanceCluster;
								find_cons_cluster_whole_cluster_distance_calculation(
									it_,  _clustersFound, cons_cluster_number_, conflict_cluster_set_, 
									uncertain_set_,  LP_Trans_Covar_Map, 
									VertexInf, 1, chiStatis, futileBit, 
									transformDistanceCluster, consCluster_serialNum);		
								cout<<"**********************111111111111111********************"<<endl;
								cout<<"trnasform distance for each element in cluster "<<toTestClusterID<<endl;
								for(int dd843= 0; dd843 < transformDistanceCluster.size(); dd843++)
								{
									cout<<transformDistanceCluster[dd843]<<" ";
								}
								cout<<endl;
							}
						}
						else
						{
							cout<<"to test cluster element number is bigger than one"<<endl;
							std::vector<double> transformDistanceCluster;
							find_cons_cluster_whole_cluster_distance_calculation(
								it,  _clustersFound, cons_cluster_number_, conflict_cluster_set_, 
								uncertain_set_,  LP_Trans_Covar_Map, 
								VertexInf, 0, chiStatis, futileBit, 
								transformDistanceCluster, toTestClusterID);	
							cout<<"*********************22222222222222222222*******************"<<endl;
							cout<<"transform distance for each element in cluster "<<toTestClusterID<<endl;
							for(int dd843= 0; dd843 < transformDistanceCluster.size(); dd843++)
							{
								cout<<transformDistanceCluster[dd843]<<" ";
							}
							cout<<endl;	

							//calculate the
							int first_local  = _clustersFound[toTestClusterID].positionserial.back()[0];
							int second_local = _clustersFound[toTestClusterID].positionserial.back()[3];		
							IntPairSet::const_iterator it_ = loops.begin(), lend_ = loops.end();

							for(;it_!=lend_;it_++)
							{
								if(((*it_).first == first_local and (*it_).second == second_local) or 
									((*it_).first == second_local and (*it_).second == first_local))
									break;
							}
							if(it_==lend_)
							{
					    		printf("can not find the loop in cluster 29. This error is in %s on line %d\n",  __FILE__, __LINE__);
					    		exit(0);
							}

							transformDistanceCluster.clear();
							find_cons_cluster_whole_cluster_distance_calculation(
								it_,  _clustersFound, cons_cluster_number_, conflict_cluster_set_, 
								uncertain_set_,  LP_Trans_Covar_Map, 
								VertexInf, 0, chiStatis, futileBit, 
								transformDistanceCluster, toTestClusterID);		
							cout<<"***************2222222*********11111111***********"<<endl;
							cout<<"trnasform distance for each element in cluster "<<toTestClusterID<<endl;
							for(int dd843= 0; dd843 < transformDistanceCluster.size(); dd843++)
							{
								cout<<transformDistanceCluster[dd843]<<" ";
							}
							cout<<endl;

							//transfrom distance between the last element of consistent cluster and other elements in the same clsuter 
							first_local  = _clustersFound[consCluster_serialNum].positionserial.back()[0];
							second_local = _clustersFound[consCluster_serialNum].positionserial.back()[3];		
							it_ = loops.begin(), lend_ = loops.end();

							for(;it_!=lend_;it_++)
							{
								if(((*it_).first == first_local and (*it_).second == second_local) or 
									((*it_).first == second_local and (*it_).second == first_local))
									break;
							}
							if(it_==lend_)
							{
					    		printf("can not find the loop in cluster 29. This error is in %s on line %d\n",  __FILE__, __LINE__);
					    		exit(0);
							}
							transformDistanceCluster.clear();
							find_cons_cluster_whole_cluster_distance_calculation(
								it_,  _clustersFound, cons_cluster_number_, conflict_cluster_set_, 
								uncertain_set_,  LP_Trans_Covar_Map, 
								VertexInf, 0, chiStatis, futileBit, 
								transformDistanceCluster, consCluster_serialNum);		
							cout<<"***************2222222*********11111111****111111111*******"<<endl;
							cout<<"trnasform distance in final consistent cluster, the base loop is back() "<<toTestClusterID<<endl;
							for(int dd843= 0; dd843 < transformDistanceCluster.size(); dd843++)
							{
								cout<<transformDistanceCluster[dd843]<<" ";
							}
							cout<<endl;
							//set each loop in final conistent cluster as base, then calculate the tranfrom distance to the both exitst cluster
							for(int eachCo = 0; eachCo < _clustersFound[consCluster_serialNum].positionserial.size(); eachCo++)
							{
								first_local  = _clustersFound[consCluster_serialNum].positionserial[eachCo][0];
								second_local = _clustersFound[consCluster_serialNum].positionserial[eachCo][3];		
								IntPairSet::const_iterator it_ = loops.begin(), lend_ = loops.end();

								for(;it_!=lend_;it_++)
								{
									if(((*it_).first == first_local and (*it_).second == second_local) or 
										((*it_).first == second_local and (*it_).second == first_local))
										break;
								}
								if(it_==lend_)
								{
									cout<<"start "<<first_local<<"  end "<<second_local<<endl;
									cout<<"(*it_).first  "<<(*it_).first<<"   (*it_).second"<<(*it_).first<<endl;
						    		printf("can not find the loop in cluster 29. This error is in %s on line %d\n",  __FILE__, __LINE__);
						    		exit(0);
								}

								transformDistanceCluster.clear();
								find_cons_cluster_whole_cluster_distance_calculation(
									it_,  _clustersFound, cons_cluster_number_, conflict_cluster_set_, 
									uncertain_set_,  LP_Trans_Covar_Map, 
									VertexInf, 0, chiStatis, futileBit, 
									transformDistanceCluster, toTestClusterID);		
								cout<<"***************2222222*******  eachCo  "<<eachCo<<"  ***********"<<endl;
								for(int dd843= 0; dd843 < transformDistanceCluster.size(); dd843++)
								{
									cout<<transformDistanceCluster[dd843]<<" ";
								}
								cout<<endl;
							}
							cout<<"***** end ******2222222**** end *****22222222****** end *****"<<endl;
						}
						if(exitBit == 1)
						{
							// sleep(2);
							// exit(0);
						}
					}
				}
				//assert(split.size() == splitDis.size());
				if(split.size() != 0)
				{
					std::pair<int, int> changeLoop;
					changeLoop.first  = _clustersFound[split[0].first].positionserial[split[0].second[0]][0];
					changeLoop.second = _clustersFound[split[0].first].positionserial[split[0].second[0]][3];

					IntPairSet::const_iterator it_ = loops.begin(), lend_ = loops.end();

					for(;it_!=lend_;it_++)
					{
						if(((*it_).first == changeLoop.first and (*it_).second == changeLoop.second) or 
							((*it_).first == changeLoop.second and (*it_).second == changeLoop.first))
							break;
					}
					std::vector<int> minDISvalid;
					find_cons_cluster_second_re( it_,  _clustersFound,  cons_cluster_number2, conflict_cluster_set,  
						uncertain_set, LP_Trans_Covar_Map, VertexInf, 1 , chiStatis,  futileBit, split2, splitDis2, minDISvalid,
						conflictClusterDistanceMap);
					if(split2.size() > 0)//if find new clster need to be splited
					{
						// splitSET, split2SET, spCONset, sp2CONse
						// splitSET.insert(split)
						// std::cout << "find new split info." << std::endl;
						// std::cout<<"split2 size: "<<split2.size()<<endl;
						// std::cout<<"split  size: "<<split.size()<<endl;
						// // assert(0);
						// if()
						// if(split2[0].first != split[0].first )
						// {
						// 	cout<<"split2 correspoding to cluster "<<split2[0].first<<endl;
						// 	std::cout<<"the cluster to previous split is "<<(split[0].first)<<std::endl;
						// 	cout<<"the minmun dis in split2 is "<<minDISvalid[0]<<endl;
						// 	if(minDISvalid[0] < 1)
						// 	{
						// 		cout<<"the minimun dis in split2 is less than 1, maybe should be considerded carefully"<<endl;
						// 		std::cin.get();	
						// 	}
						// }
						// else
						// 	cout<<"the cluster to split is the same one"<<endl;
						bool find = 0;
						for(int split2Num = 0; split2Num < split2.size(); split2Num++)
						{
							find = 0;
							for(int split1Num = 0; split1Num < split.size(); split1Num++)
							{
								if(split2[split2Num].first == split[split1Num].first)
								{
									find = 1;
								}
							}
							if(!find)
							{
								// minDistance2 = ;
								cout<<"find new cluster that has to be splited, whose minimun distance is "<<minDISvalid[split2Num]<<endl;
								vectorNewSplit.push_back(minDISvalid[split2Num]);
								if(minDISvalid[split2Num] < 1)
								{
									vectorNewSplitCloser.push_back(minDISvalid[split2Num]);
									// std::cin.get();	
								}
							}
						}
					}	
					//if find new consistent info, find out the relaitonship with the split	
					int consistentClusterToHoldSplit = -1;	
					double disToPickBestCluster = 100.1;			
					if(cons_cluster_number2.size() > 0)
					{
						for(int split2Num = 0; split2Num < cons_cluster_number2.size(); split2Num++)
						{
							bool retURN = 0;
							for(int split1Num = 0; split1Num < split.size(); split1Num++)
							{
								if(cons_cluster_number2[split2Num] == split[split1Num].first)
								{
									cout<<"this consistent cluster "<<cons_cluster_number2[split2Num]<<" equal to first time split cluster "
										<<split[split1Num].first<<endl;
									printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
									retURN = 1;
									// std::cin.get();
									break;	
								}


							}
							if(!retURN )
							{
								if(chiStatis[split2Num] < disToPickBestCluster)
								{
									consistentClusterToHoldSplit = cons_cluster_number2[split2Num];
									disToPickBestCluster = chiStatis[split2Num];
								}
							}
						}
						// cout<<"consistent cluster number in second procedure is "<<cons_cluster_number2.size()<<endl;
						// cout<<"consistent cluster id is "<<cons_cluster_number2[0]<<endl;
						// if(cons_cluster_number2.size() > 1)
						// {
						// 	cout<<"find more than one cosistent cluster when determine the relationship of split"<<endl;
						// 	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						// 	std::cin.get();	
						// }
						// if(cons_cluster_number2.size() == 1 and cons_cluster_number2[0] == split[0].first)
						// {
						// 	if(split.size() > 0)
						// 	{
						// 		cout<<"the minimun dis in split2 is less than 1, maybe should be considerded carefully"<<endl;
						// 		printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						// 		std::cin.get();	
						// 	}
						// }
						// else
						// {
						// 	std::cout << "find consistent clusters to the split cluster." << std::endl;
						// 	for(int toPrint = 0; toPrint < cons_cluster_number2.size(); toPrint++)
						// 		std::cout<<cons_cluster_number2[toPrint]<<" ";
						// 	std::cout<<endl;

						// 	for(int toPrint = 0; toPrint < split2.size(); toPrint++)
						// 	{
						// 		cout<<"* * the cluster need to be splited is "<<split2[toPrint].first<<endl;
						// 		for(int toPrintSecond = 0; toPrintSecond < split2[toPrint].second.size(); toPrintSecond++)
						// 		{
						// 			cout<<split2[toPrint].second[toPrintSecond]<<" ";
						// 		}
						// 		cout<<endl;
						// 	}
						// 	// if split info in two process is different
						// 	std::cout<<"the cluster to previous split is "<<(split[0].first)<<std::endl;
						// 	std::cin.get();	
						// } updateSplit_original
					}

					updateSplit(split, splitDis, _clustersFound, cons_cluster_number2, 
						conflict_cluster_set, all_conflict_cluster, all_uncertain_cluster, consistentClusterToHoldSplit);

					//add element to the cluster pair adn distance mapp
					if(conflict_cluster_set.size() + conflictClusterDistanceMap.size() > 0)
					{
						if(conflict_cluster_set.size() != conflictClusterDistanceMap.size())
						{
							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);		
						}

						int toHold = 0;
						if(consistentClusterToHoldSplit != -1)
							toHold = consistentClusterToHoldSplit;
						else
							toHold = _clustersFound.size() - 1;


						for(int addMap = 0; addMap < conflict_cluster_set.size(); addMap++)
						{
							conflictClusterPairDistanceMap[std::pair<int,int> (conflict_cluster_set[addMap], toHold)] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];
							conflictClusterPairDistanceMap[std::pair<int,int> (toHold, conflict_cluster_set[addMap])] =
							 	conflictClusterDistanceMap[conflict_cluster_set[addMap]];  
						}
					}

				}



				// if((start == 1590  and end == 25) or ((end == 1590  and  start == 25)))
				// {
				// 		std::cout << "start == 1590  and end == 25, Press \'Return\' to end." << std::endl;
				// 		cout<<"consCluster_serialNum: "<<consCluster_serialNum<<endl;
				// 		std::cin.get();					
				// }

				// if(consCluster_serialNum == 13)
				// 	{
				// 		std::cout << "consCluster_serialNum == 13, Press \'Return\' to end." << std::endl;
				// 		std::cin.get();
				// 	}

			}
			// if(consCluster_serialNum == 13)
			// // if((*it).first ==9016 and (*it).second ==3876)
			// {
			// 	cout<<"get the loop's id when add to cluster "<<loopToClusterIDMap[*it]<<endl;
			// 	std::cin.get();
			// 	// exit(0);
				
			// 	// int itmc = find_ele.find(all_conflict_cluster,  consCluster_serialNum);
			// 	// if(itmc == -1)
			// 	// {
			// 	// 	cerr<<"the size of all conflict info: "<<all_conflict_cluster.size()<<endl;
			// 	// 	cerr<<"the elemets in all conflict info: ";
			// 	// 	for(int local_de = 0; local_de < all_conflict_cluster.size(); local_de++)
			// 	// 		cerr<<all_conflict_cluster[local_de].first<<" ";
			// 	// 	cerr<<endl;
			// 	// 	cerr<<"failed to find the conflict info"<<endl;
		 //  //   		// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		 //  //   		// exit(0);
			// 	// }
			// }
			// if(_clustersFound.size() > 85)
			// {
			// 	cout<<"loopToClusterIDMap[std::pair<int,int> (9016, 3876)]: "<<loopToClusterIDMap[std::pair<int,int> (9016, 3876)]<<endl;
			// 	if(loopToClusterIDMap[std::pair<int,int> (9016, 3876)] > 200)
			// 		cout<<"try to find why cluster id abnormal "<<loopToClusterIDMap[std::pair<int,int> (9016, 3876)]<<endl;
			// 	std::cin.get();	
			// }

			// 	if(itmc != -1)
			// 	{
			// 		cout<<"cluster "<<consCluster_serialNum<<"'s conflict cluster are: ";
			// 		for(int pp = 0; pp < all_conflict_cluster[itmc].second.size(); pp++)
			// 		{
			// 			cout<<" "<<all_conflict_cluster[itmc].second[pp];
			// 		}
			// 		cout<<endl;					
			// 	}


			// 	itmc = find_ele.find(all_uncertain_cluster,  consCluster_serialNum);
			// 	if(itmc == -1)
			// 	{
			// 		cout<<"failed to find the uncertain info"<<endl;
		 //    		printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		 //    		// exit(0);				
			// 	}
			// 	else
			// 	{
			// 		cout<<"cluster "<<consCluster_serialNum<<"'s uncertain cluster are: ";
			// 		for(int pp = 0; pp < all_uncertain_cluster[itmc].second.size(); pp++)
			// 		{
			// 			cout<<" "<<all_uncertain_cluster[itmc].second[pp];
			// 			if(all_uncertain_cluster[itmc].second[pp] == -1)
			// 			{
			// 	    		printf("find 210. This error is in %s on line %d\n",  __FILE__, __LINE__);
			// 	    		exit(0);		
			// 			}
			// 		}
			// 		cout<<endl;
			// 	}

			// 	std::cout << "consCluster_serialNum = "<<consCluster_serialNum<<" , Press \'Return\' to continue." << std::endl;
			// 	// std::cin.get();
			// }

		}

		std::array<double,6> ty={1,1,1,1,1,1};

		// cout<<"get the loop's id when add to cluster "<<loopToClusterIDMap[std::pair<int,int> (9016, 3876)]<<endl;
		// std::cin.get();

		// fileStream.open("clusterFile.g2o",ios::trunc);
		ofstream fileStreamr; 
		fileStreamr.open(nameofclusterfile,ios::trunc);
	
		// // cout<<"clusters:"<<endl;
		// fileStreamr.open("cluster-all.txt",ios::trunc);
		std::pair<g2o::SE2, Matrix3d> tSave;
		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			fileStreamr<<i<<"\n";
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[i].positionserial.begin(), 
				lendVertex = _clustersFound[i].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				tSave = LP_Trans_Covar_Map[std::pair<int,int> (ty[0], ty[3])];

				fileStreamr<<"EDGE_SE2 "<<ty[0]<<" "<<ty[3]<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
					<<" "<<1.0/tSave.second(0,0)<<" "<<0<<" "<<0<<" "<<1.0/tSave.second(1,1)<<" "<<0<<" "<<1.0/tSave.second(2,2)<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();

		//display the conflict infor of each cluster
		cout<<"confilict info: "<<endl;
		for(int itmc = 0; itmc < all_conflict_cluster.size(); itmc++)
		{
			cout<<"cluster "<<all_conflict_cluster[itmc].first<<"'s conflict cluster are: ";
			for(int j = 0; j < all_conflict_cluster[itmc].second.size(); j++)
			{
				cout<<" "<<all_conflict_cluster[itmc].second[j];
			}
			cout<<endl;
		}

		cout<<" "<<endl;
		cout<<"uncertain info: "<<endl;

		for(int itmc = 0; itmc < all_uncertain_cluster.size(); itmc++)
		{
			cout<<"cluster "<<all_uncertain_cluster[itmc].first<<"'s uncertain cluster are: ";
			for(int j = 0; j < all_uncertain_cluster[itmc].second.size(); j++)
			{
				cout<<" "<<all_uncertain_cluster[itmc].second[j];
			}
			cout<<endl;
		}

		//check the info in conflich and that in uncertian set to see if there are some elements exit in both sets
		cout<<"check uncertain and confilict to see if there is a element exits in both sets: "<<endl;
		int findNoAccordingInUncertianSet;
		int findsecond;
		for(int itmc = 0; itmc < all_conflict_cluster.size(); itmc++)
		{
			cout<<" "<<endl;
			cout<<"current cluster  in confict set is "<<all_conflict_cluster[itmc].first<<endl;

			findNoAccordingInUncertianSet = find_ele.find(all_uncertain_cluster,  all_conflict_cluster[itmc].first);

			if(findNoAccordingInUncertianSet == -1)
			{
				cout<<"no according uncertian info to that cluster"<<endl;
			}
			else
			{
				cout<<"find according uncertian info to that cluster"<<all_conflict_cluster[itmc].first<<endl;
				
				for(int j = 0; j < all_conflict_cluster[itmc].second.size(); j++)
				{
					findsecond = find_ele.find(all_uncertain_cluster[findNoAccordingInUncertianSet].second, all_conflict_cluster[itmc].second[j]);
					if(findsecond != -1)
					{
						cout<<"cluster "<<all_conflict_cluster[itmc].first<<"'s conflict cluster are: ";
						for(int pp = 0; pp < all_conflict_cluster[itmc].second.size(); pp++)
						{
							cout<<" "<<all_conflict_cluster[itmc].second[pp];
						}
						cout<<endl;
						cout<<"cluster "<<all_uncertain_cluster[findNoAccordingInUncertianSet].first<<"'s uncertain cluster are: ";
						for(int pp = 0; pp < all_uncertain_cluster[findNoAccordingInUncertianSet].second.size(); pp++)
						{
							cout<<" "<<all_uncertain_cluster[findNoAccordingInUncertianSet].second[pp];
						}
						cout<<endl;						
						// cout<<"cluster "<<all_conflict_cluster[itmc].second[j]<<" exits in both confilict and uncertain set"<<endl;
				  //   	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				  //   	exit(0);
					}
				}
			}
		} 

		for(int itmc = 0; itmc < all_uncertain_cluster.size(); itmc++)
		{
			cout<<" "<<endl;
			cout<<"current cluster  in uncertain set is "<<all_uncertain_cluster[itmc].first<<endl;

			findNoAccordingInUncertianSet = find_ele.find(all_conflict_cluster ,  all_uncertain_cluster[itmc].first);

			if(findNoAccordingInUncertianSet == -1)
			{
				cout<<"no according conflict info to that cluster"<<endl;
			}
			else
			{
				cout<<"find according conflict info to that cluster"<<all_uncertain_cluster[itmc].first<<endl;
				// for(int j = 0; j < all_uncertain_cluster[itmc].second.size(); j++)
				// {
				// 	findsecond = find_ele.find(all_conflict_cluster[findNoAccordingInUncertianSet].second, all_uncertain_cluster[itmc].second[j]);
				// 	if(findsecond != -1)
				// 	{
				// 		cout<<"cluster "<<all_uncertain_cluster[itmc].second[j]<<" exits in both confilict and uncertain set"<<endl;
				//     	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				//     	exit(0);
				// 	}
				// }
			}

		} 

		// std::cout << "after clustering, Press \'Return\' to end." << std::endl;
		// std::cin.get();

#if 0
		if(0)
		{
			std::cout<<" \% Clusters formed "<<_clustersFound.size()<<std::endl;
			std::cout<<"limits = [ "<<std::endl;
			for(size_t i=0 ; i< _clustersFound.size() ; i++)
			{
				std::cout<<i<<" -> sz "<<_clustersFound[i].size<<" :: ";
				std::cout<<" "<<_clustersFound[i].startLow<<" "<<_clustersFound[i].startHigh<<" ";
				std::cout<<" "<<_clustersFound[i].endLow<<" "<<_clustersFound[i].endHigh<<std::endl;;

			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;


			std::cout<<"membership =[ ";
			for(size_t i=0; i<membership.size();i++)
			{
				std::cout<<membership[i]<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;
		}
#endif

	}

// 	/* collect vertex ID and position and angle in vertexInfo
	
// 	 */

	int collect_vertexAndLP(const char* filename, std::map<std::pair<int, int>, std::array<double,9> > &LC_Inf,
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map)
	{

		ifstream fileStream;  

	    string tmp,temp1;  
	    std::array<double,4> verT={0, 0, 0, 0};
	    std::array<double,11> odoedge_element;
	    std::array<double,11> savemid;
	    std::array<double,9> lcedge_element;
	    std::pair<int, int>  lc_vertex_pair;
	    std::pair<g2o::SE2, Matrix3d> ele_lp_trans_covar;
	    // char* seg;
	    int count = 0,dddlndgsfdgj=0;// 行数计数器  
	    int position = 0, position2 = 0;  
	    double nul;

	    // fileStream.open("B25b_0.500.g2o",ios::in);//ios::in 表示以只读的方式读取文件
	    fileStream.open(filename,ios::in);  
	    if(fileStream.fail())//文件打开失败:返回0  
	    { 
	    	cout<<"open file failed"<<endl; 
	    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
	    	exit(0);
	        return 0;  
	    }  
	    else//文件存在  
	    {  
		    while(getline(fileStream,tmp,'\n'))//读取一行  
		    {  
		    	position = tmp.find("VERTEX_SE2");
		    	position2 = tmp.find("EDGE_SE2");
		    	// cout<<tmp<<endl; 	
		    	// exit(0);
		
		    	istringstream stringin(tmp);
		    	if (tmp.size() > 0 )
		    	{
					if(position!=string::npos)  
			    	{			    		
			    		for (dddlndgsfdgj=0;dddlndgsfdgj<5;dddlndgsfdgj++)
			    		{
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1:
			    				case 2:
			    				case 3:
			    				case 4:
			    					stringin >> verT[dddlndgsfdgj-1];
			    					break;
			    				default:
			    					break;
			    			}
			    		}
			    	
			    		VertexInf.push_back(verT);
			    		// cout<<verT[0]<<" "<<verT[1]<<" "<<verT[2]<<" "<<verT[3]<<endl;
			    		// cout<<VertexInf.size()<<endl;
			    		// cout<<VertexInf[0][0]<<" "<<VertexInf[0][1]<<" "<<VertexInf[0][2]<<" "<<VertexInf[0][3]<<endl;
			    		// exit(0);
			    	}
			    	else if(position2 !=string::npos)
			    	{
			    		bool odobit = 0;

						for (dddlndgsfdgj = 0; dddlndgsfdgj < 12; dddlndgsfdgj++)
			    		{	
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1://start vertex numver
			    					stringin >> lc_vertex_pair.first;
			    					break;
			    				case 2://end vertex number
			    					stringin >> lc_vertex_pair.second;
			    					if(std::abs(lc_vertex_pair.first - lc_vertex_pair.second) == 1)
			    						odobit = 1;
			    					break;
			    				case 3://x of the edge
			    					if(odobit)//if its a odometry edge put in vector[8]
			    					{
			    						odoedge_element[0] = lc_vertex_pair.first;
			    						odoedge_element[1] = lc_vertex_pair.second;
			    						stringin >> odoedge_element[2];
			    					}
			    					else
			    						stringin >> lcedge_element[0];
			    					break;
			    				case 4: //y of the edge
			    				case 5: //andlge of the adge
			    				case 6: //information matrix (0,0)	
			    				case 7: //information matrix (0,1)	
			    				case 8: //information matrix (0,2)	
			    				case 10://information matrix (1,1)	
			    				case 9: //information matrix (1,2)
			    				case 11://information matrix (2,2)
			    					if(odobit)//if its a odometry edge put in vector[8]
			    						stringin >> odoedge_element[dddlndgsfdgj-1];
			    					else
			    						stringin >> lcedge_element[dddlndgsfdgj-3];
			    					break;	
			    				default:
			    				{
			    					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
									exit(0);
			    					break;
			    				}
			    			}
			    		}
			    		
			    		if(odobit)
			    			OdoInf.push_back(odoedge_element);
			    		else
			    		{
			    			count++; 
			    			LC_Inf[lc_vertex_pair] = lcedge_element;
				    		// cout<<lc_vertex_pair.first<<" "<<lc_vertex_pair.second<<" "<<lcedge_element[0]<<" "
				    		// 	<<lcedge_element[1]<<" "<<lcedge_element[2]<<" "<<lcedge_element[3]<<" "<<lcedge_element[4]<<" "
			    			// 	<<lcedge_element[5]<<endl;	

							Matrix3d m1 = Matrix3d::Identity();
							g2o::Vector3 mid_vector3; 
							g2o::SE2 edge1;

							m1(0,0) = lcedge_element[3]; 
							m1(0,1) = lcedge_element[4]; 
							m1(1,0) = lcedge_element[4]; 

							m1(0,2) = lcedge_element[5];
							m1(2,0) = lcedge_element[5];

							m1(1,1 )= lcedge_element[6]; 

							m1(1,2) = lcedge_element[7]; 
							m1(2,1) = lcedge_element[7]; 

							m1(2,2) = lcedge_element[8];


							Matrix3d mx;
							mx = m1.inverse();

							mid_vector3[0] = lcedge_element[0];
							mid_vector3[1] = lcedge_element[1];
							mid_vector3[2] = lcedge_element[2];
							edge1.fromVector(mid_vector3);

							ele_lp_trans_covar.first  = edge1;
							ele_lp_trans_covar.second = mx;

							LP_Trans_Covar_Map[lc_vertex_pair] = ele_lp_trans_covar;
							// if(lc_vertex_pair.first == 1575)
							// {
							// 	cout<<"lc_vertex_pair.first: "<<lc_vertex_pair.first<<endl;

							// 	exit(0);					
							// }
							// else if(lc_vertex_pair.first > 1575)
							// {
							// 	cout<<"lc_vertex_pair.first: "<<lc_vertex_pair.first<<endl;
							// 	exit(0);					
							// }

			    			if(lc_vertex_pair.first <lc_vertex_pair.second){
			    				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
								exit(0);
			    			}
			    		}

			    	}
		    	} 
		    } 

		    // cout<<"loop num: "<<LP_Trans_Covar_Map.size()<<endl;
		    // exit(0);

		    if(OdoInf.size() != VertexInf.size()-1)
		    {
		    	cout<<"odo size: "<<OdoInf.size()<<" vertex size: "<<VertexInf.size()<<endl;
		    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    for(int i=0;i<OdoInf.size();i++)
		    {
		    	if(std::abs(OdoInf[i][0] - OdoInf[i][1]) != 1 )
		    	{
		    		cout<<"OdoInf[i][0]: "<<OdoInf[i][0]<<" OdoInf[i][1]: "<<OdoInf[i][1]<<endl;
		    		cout<<"i: "<<i<<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
		    	}
		    	else if(OdoInf[i][0] != i)
		    	{
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
		    	}
		    	// cout<<OdoInf[i][0]<<" "<<OdoInf[i][1]<<" "<<OdoInf[i][2]<<" "<<OdoInf[i][3]<<" "<<OdoInf[i][4]<<" "
		    	// 	<<OdoInf[i][5]<<" "<<OdoInf[i][6]<<" "<<OdoInf[i][7]<<endl;
		    }

		    printf("loop count:%d \n",count);
		    printf("sum:%d \n",int(VertexInf.size()+LC_Inf.size()+OdoInf.size()));//+LC_Inf.size() +OdoInf.size()
		    cout<<"VertexInf.size():"<<VertexInf.size()<<endl;
		    fileStream.close();  
		    // exit(0);

		    // if(filename.find() != )
		    // 	tmp.find("VERTEX_SE2")
		    cout<<filename<<endl;
		    string calLineOfTLC, TLCFileName;  
		    calLineOfTLC = filename;
		    cout<<calLineOfTLC<<endl;

		    int ret3 = calLineOfTLC.find("_alpha");
		    cout<<calLineOfTLC<<endl;
		    cout<<calLineOfTLC.find("_alpha")<<endl;
		    if(ret3 != -1)
		    {
		    	NumOfRealLoops;//number of the true loop closures , namely the line numbers in the txt file
		    	int pos = calLineOfTLC.find(".g2o");
		    	if(pos == -1)
		    	{
		    		cout<<"can not find \' .g2o\' in input file name ,so cant't calculate the number of true loops"<<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
		    	}
		    	TLCFileName = calLineOfTLC;//delete the last four chars and then add the part of GT LC.txt
		    	TLCFileName.replace(pos,4,"_GT_LC.txt");
 				
 				fileStream.open(TLCFileName,ios::in);  
			    if(fileStream.fail())//文件打开失败:返回0  
			    {  
			    	cout<<"GT Loop txt file exists but can not be open !"<<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			    	exit(0);
			    }  
			    else//文件存在  
			    {  
			    	NumOfRealLoops = 0;
			    	std::pair<int, int> ele_GTL;
				    while(getline(fileStream,tmp,'\n'))//读取一行  
				    { 
				    	istringstream stringin(tmp);
				    	if (tmp.size() > 0 )
				    	{
				    		// int i = tmp.size();
				    		// cout<<"tmp size: "<<i<<endl;
				    		// for(int j =0; j <i; j++)
				    		// {
				    		// 	cout<<"tmp["<<j<<"] = "<<tmp[j]<<endl;
				    		// }
				    		
					    		for (int i=0; i < 2;i++)
					    		{
					    			// cout<<"tmp["<<i<<"] = "<<tmp[i]<<endl;
					    			switch(i)
					    			{
				                    	case 0:
					    					stringin >> ele_GTL.first;
					    					// cout<<temp1<<endl;
					    					break;
					    				case 1:
					    					stringin >> ele_GTL.second;
					    					break;
					    			}
					    		}
					    }
					    // exit(0);
				    	set4GTL.insert(ele_GTL);
						NumOfRealLoops = NumOfRealLoops+1;
				    } 		    	
				}
				fileStream.close(); 
			} 
			if(NumOfRealLoops != set4GTL.size())
			{
				if(NumOfRealLoops == -1)
					cout<<"no ground truth LC file available"<<endl;
				else
				{
			    	cout<<"the number of true loops does not consistent!"<<endl;
			    	cout<<"NumOfRealLoops: "<< NumOfRealLoops<<" set4GTL.size(): "<<set4GTL.size() <<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			    	exit(0);	
				}
			}
		    // calLineOfTLC
		    cout<<ret3<<endl; 
		    cout<<"final name: "<<TLCFileName<<" has "<< NumOfRealLoops<<" true loops."<<endl;
		    // exit(0);
		    return count ;  
		} 
	}

	void prepare_accelerateBySynthetic(std::vector<std::array<double,11>> & OdoInf, 
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo100_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo1000_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10000_trans_varMatrix,
		std::vector<double > & syntheticOdo10_dis, std::vector<double > & syntheticOdo100_dis,
		std::vector<double > & syntheticOdo1000_dis, std::vector<double > & syntheticOdo10000_dis)
	{
		int start, end;
		std::pair<g2o::SE2, Matrix3d> result_synthetic_odo;
		std::array<double, 5>  length;
		bool  fultileBit;
		for(int i = 0; i*10+10 <= OdoInf.size(); i++ )
		{
			start = i*10;
			end   = i*10+10;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo10_dis.push_back(length[0]);
		}
		for(int i = 0; i*100+100 <= OdoInf.size(); i++ )
		{
			start = i*100;
			end   = i*100+100;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo100_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo100_dis.push_back(length[0]);
		}
		for(int i = 0; i*1000+1000 <= OdoInf.size(); i++ )
		{
			start = i*1000;
			end   = i*1000+1000;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo1000_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo1000_dis.push_back(length[0]);
		}
		for(int i = 0; i*10000+10000 <= OdoInf.size(); i++ )
		{
			start = i*10000;
			end   = i*10000+10000;
			synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first,
				result_synthetic_odo.second, length, fultileBit);
			syntheticOdo10000_trans_varMatrix.push_back(result_synthetic_odo);
			syntheticOdo10000_dis.push_back(length[0]);
		}
	}

	void second_prepare_accelerated_synthetic_odo(std::vector<std::array<double,11>> & OdoInf, 
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdoAccumulate_trans_varMatrix,
		std::vector<double> & syntheticOdoAccumulate_dis)
	{
		int i = 0;
		std::array<double, 4> length_ele;
		std::pair<g2o::SE2, Matrix3d>  result_synthetic_odo, result_synthetic_odo2;
		std::array<double, 5>  length;
		bool fultileBit;
		int odoSize = OdoInf.size();

		while(i <= OdoInf.size())
		{
			
			accelerated_synthetic_odo(0,  i, OdoInf, 
			result_synthetic_odo,  length, syntheticOdo10_trans_varMatrix,
			syntheticOdo100_trans_varMatrix,   syntheticOdo1000_trans_varMatrix,
			syntheticOdo10000_trans_varMatrix, syntheticOdo10_dis, syntheticOdo100_dis,
			syntheticOdo1000_dis, syntheticOdo10000_dis, fultileBit);

			syntheticOdoAccumulate_dis.push_back(length[0]);
			syntheticOdoAccumulate_trans_varMatrix.push_back(result_synthetic_odo);

			synthesize_odo_edges( 0, i, OdoInf, result_synthetic_odo2.first, 
				result_synthetic_odo2.second, length, fultileBit);

			judgeTwoSE2Equal(result_synthetic_odo.first, result_synthetic_odo2.first);

			if(judgeTwoMatrix3dEqual(result_synthetic_odo.second,result_synthetic_odo2.second))
			{
				cout<<"start: 0 "<<" end: "<<i<<endl;
				cout<<"cov2:"<<endl;
				cout<<result_synthetic_odo.second<<endl;
				cout<<"cov1_:"<<endl;
				cout<<result_synthetic_odo2.second<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);			
			}
			i++;
		}
		// ;
		// result_synhetic_odo2.second;

	}



	void second_accelerated_synthetic_odo(int start, int end,
		std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::array<double, 5>  & length,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdoAccumulate_trans_varMatrix,
		std::vector<double > & syntheticOdoAccumulate_dis,
		bool & fultileBit)
	{	
		g2o::SE2 midTrans, nextTrans, SE; 
		Matrix3d nextVarMatrix, midVarMatrix, j1, j2, cov;		

		g2o::SE2 toFindOutWhyAcceleratedNotEqualSlow_SE2;
		Matrix3d toFindOutWhyAcceleratedNotEqualSlow_VAR;
		bool toFindOutWhyAcceleratedNotEqualSlow_fultileBit;
		std::array<double,5> toFindOutWhyAcceleratedNotEqualSlow_length;

		if(start == 0)
		{
			result_synthetic_odo = syntheticOdoAccumulate_trans_varMatrix[end];
			SE = syntheticOdoAccumulate_trans_varMatrix[end].first;
		}
		else
		{
			midTrans     = syntheticOdoAccumulate_trans_varMatrix[start].first.inverse();
			SE = midTrans;
			midVarMatrix = syntheticOdoAccumulate_trans_varMatrix[start].second;
		    
		    // cerr<<"vitualStart: "<<vitualStart <<endl;
			Jacobian_4_edge_propagate(midTrans, syntheticOdoAccumulate_trans_varMatrix[end].first, j1, j2);
			covariance_propagate(midVarMatrix, syntheticOdoAccumulate_trans_varMatrix[end].second, j1, j2, nextVarMatrix);
			result_synthetic_odo.first = midTrans*(syntheticOdoAccumulate_trans_varMatrix[end].first);
			result_synthetic_odo.second = nextVarMatrix;
			length[0] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
			length[2] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
			length[3] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];
			length[4] = syntheticOdoAccumulate_dis[end] - syntheticOdoAccumulate_dis[start];			
		}
		// fultileBit = 0;
		fultileBit = judgeVarianceAndDistance_Small(length);

		cout<<"midVarMatrix: "<<endl;
		cout<<midVarMatrix;
		cout<<"nextVarMatrix: "<<endl;	
		cout<<nextVarMatrix<<endl;
		
	}

	void accelerated_synthetic_odo(int start, int end, std::vector<std::array<double,11>> & OdoInf, 
		std::pair<g2o::SE2, Matrix3d> & result_synthetic_odo, std::array<double, 5>  & length,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo100_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo1000_trans_varMatrix,
		std::vector<std::pair<g2o::SE2, Matrix3d> > & syntheticOdo10000_trans_varMatrix,
		std::vector<double > & syntheticOdo10_dis, std::vector<double > & syntheticOdo100_dis,
		std::vector<double > & syntheticOdo1000_dis,  std::vector<double > & syntheticOdo10000_dis,
		bool & fultileBit)
	{	
		g2o::SE2 midTrans, nextTrans; 
		Matrix3d nextVarMatrix, midVarMatrix, j1, j2,j3, cov;		
		int a ,b ,c;

		if (start > 99999 or end > 99999)
		{
			cout<<"the number of vectex is too big"<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(start > OdoInf.size() or end > OdoInf.size())
		{
			cout<<"start: "<<start<<" end:"<<end <<" OdoInf.size:"<<OdoInf.size() <<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		int startWan = 0, startQian= 0, startBai = 0, startShi = 0, startGe  = 0;
		int endWan   = 0,  endQian = 0,  endBai  = 0,  endShi  = 0,  endGe   = 0;

		if(start > 9999)
		{
			startWan =  start / 10000;          startQian= (start % 10000) /1000;
			startBai = (start % 1000) /100;     startShi = (start % 100) /10;       startGe  =  start % 10;
		}
		else if(start > 999)
		{
			startQian= (start % 10000) /1000;   startBai = (start % 1000) /100;
			startShi = (start % 100) /10;       startGe  =  start % 10;
		}
		else if(start > 99)
		{
			startBai = (start % 1000) /100;     startShi = (start % 100) /10;      startGe  =  start % 10;
		}
		else if(start > 9)
		{
			startShi =  (start % 100) /10;      startGe  =  start % 10;
		}
		else
			startGe  =  start;			


		if(end > 9999)
		{
			endWan =  end / 10000;         endQian= (end % 10000) /1000;
			endBai = (end % 1000) /100;    endShi = (end % 100) /10;       endGe  =  end % 10;
		}
		else if(end > 999)
		{
			endQian= (end % 10000) /1000;  endBai = (end % 1000) /100;
			endShi = (end % 100) /10;      endGe  =  end % 10;
		}
		else if(end > 99)
		{
			endBai = (end % 1000) /100;    endShi = (end % 100) /10;       endGe  =  end % 10;
		}
		else if(end > 9)
		{
			endShi =  (end % 100) /10;     endGe  =  end % 10;
		}
		else
			endGe  =  end;		

		int vitualStart = start , virtualEnd =  end;	
		g2o::SE2 toFindOutWhyAcceleratedNotEqualSlow_SE2;
		Matrix3d toFindOutWhyAcceleratedNotEqualSlow_VAR;
		bool toFindOutWhyAcceleratedNotEqualSlow_fultileBit;
		std::array<double,5> toFindOutWhyAcceleratedNotEqualSlow_length;

		// cout<<"start: "<<start<<endl;
		//ge
		if ((start - startGe +20) > end) 
		    synthesize_odo_edges( start, end, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
		//shi
		else if ((start - startGe- startShi*10 +200) > end) 
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			vitualStart = increaseToTen;
			// cerr<<"vitualStart: "<<vitualStart <<endl;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);	
			int decreaseToTen;
			decreaseToTen = (end - endGe)/10 -1;
			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

			//from start to the nearest  integel multiple of ten
		    for(int tenSerial = (increaseToTen) / 10; tenSerial <= decreaseToTen; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 10;  	
		    	// cerr<<"vitualStart: "<<vitualStart <<endl;
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1;   
		    	// cout<<"vitualStart: "<<vitualStart <<endl;
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length);
		}
		//bai wei 
		else if ((start - startGe- startShi*10 - startBai*100+2000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			//bai ge
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, 
				result_synthetic_odo.second, length, fultileBit);

			vitualStart = increaseToTen; 
			// cout<<"vitualStart: "<<vitualStart <<endl;	

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

			//bai shi from start to the nearest  integel multiple of ten
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to bai
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((end -endGe - endShi*10) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;  	

				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to shi
		    a = (end - endGe - endShi*10) /10;
		    b = ((end - endGe) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;  	
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }

		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1; 
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    fultileBit = judgeVarianceAndDistance_Small(length);
		}
		//thousand 
		else if ((start - startGe- startShi*10 - startBai*100 - startQian*1000+20000) > end)
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			vitualStart = increaseToTen;
			// cout<<"start: "<<start<<endl;
			// cout<<"ge shi bai qian: "<<startGe<<" "<<startShi<<" "<<startBai<<" "<<startQian<<endl;
			// cout<<"vitualStart: "<<vitualStart <<endl;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			// cout<<"result_synthetic_odo.second: "<< result_synthetic_odo.second<<endl;
			// synthesize_odo_edges( start, vitualStart, OdoInf, toFindOutWhyAcceleratedNotEqualSlow_SE2, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_VAR, toFindOutWhyAcceleratedNotEqualSlow_length, 
			// 	toFindOutWhyAcceleratedNotEqualSlow_fultileBit);

			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to qian
		    a = (start -startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"a: "<<a<<endl; 
		    	// cout<<"start: "<<start<<"startGe: "<<startGe<<"startShi: "<<startShi<<"start: "<<start<<"start: "<<start<<endl; 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);

		    }
		    //bai to shi
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;  
		    	// cout<<"vitualStart: "<<vitualStart <<endl;	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
		    	cout<<"vitualStart: "<<vitualStart <<endl;
		    	cout<<"end: "<< end<<endl;
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }

		    fultileBit = judgeVarianceAndDistance_Small(length);
		}
		//wan
		else
		{
			int increaseToTen ;
			if(startGe != 0)
				increaseToTen = start - startGe +10;
			else
				increaseToTen = start;
			synthesize_odo_edges( start, increaseToTen, OdoInf, result_synthetic_odo.first, result_synthetic_odo.second, length, fultileBit);
			// int a = (end - endGe -(start - startGe +10)) /10;
			int end_ten = (end - endGe) /10 -1;

			//should be an integer multiple of ten
			if((end - endGe -(start - startGe +10)) % 10 != 0)
			{
				cout<<"the number of vectex is too big"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			vitualStart = increaseToTen;
			// cout<<"vitualStart: "<<vitualStart <<endl;
			//shi to bai
			a = increaseToTen / 10;
			b = (start -startGe - startShi*10 +100) /10 -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 10;		
		    	// cout<<"vitualStart: "<<vitualStart <<endl;						
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to qian
		    a = (start -startGe - startShi*10 +100) /100;
		    b = ((start -startGe - startShi*10 - startBai*100 +1000 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {   
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo100_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		   
		    //qian to wan
		    a = (start - startGe - startShi*10 - startBai*100 +1000) /1000;
		    b = ((start - startGe - startShi*10 - startBai*100 - startQian*1000 + 10000) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	 	
		    	// cout<<"vitualStart: "<<vitualStart <<endl;							
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //wan to wan
		    a = (start - startGe - startShi*10 - startBai*100 - startQian*1000 +10000) /10000;
		    b = ((end - endGe - endShi*10 - endBai*100 - endQian*1000) /10000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10000;	 
		    	// cout<<"vitualStart: "<<vitualStart <<endl;								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }

		    //wan to qian
		    a = (end - endGe - endShi*10 - endBai*100 - endQian*1000) /1000;
		    b = ((end -endGe - endShi*10 - endBai*100) /1000) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 1000;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo1000_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo1000_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo1000_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //qian to bai
		    a = ((end -endGe - endShi*10 - endBai*100) /100) ;
		    b = ((end -endGe - endShi*10 ) /100) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 100;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo100_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo100_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo100_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //bai to shii
		    a = ((end -endGe - endShi*10) /10) ;
		    b = ((end -endGe ) /10) -1;
		    for(int tenSerial = a; tenSerial <= b; tenSerial++)
		    {  
		    	vitualStart = vitualStart + 10;	
		    	// cout<<"vitualStart: "<<vitualStart <<endl; 								
				Jacobian_4_edge_propagate(result_synthetic_odo.first, syntheticOdo10_trans_varMatrix[tenSerial].first, j1, j2);
				covariance_propagate(result_synthetic_odo.second, syntheticOdo10_trans_varMatrix[tenSerial].second, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * (syntheticOdo10_trans_varMatrix[tenSerial].first);
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + syntheticOdo10_dis[tenSerial];
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    //from the last  integel multiple of ten to end
		    g2o::Vector3 mid_vector3;
		    std::array<double,11> startNodeInfo;
		    for(int tenSerial = end - endGe; tenSerial < end; tenSerial++)
		    { 
		    	vitualStart = vitualStart + 1;
		    	// cout<<"vitualStart: "<<vitualStart <<endl;  	
		    	startNodeInfo = OdoInf[tenSerial];

				cov(0,0) = startNodeInfo[5]; 
				cov(0,1) = startNodeInfo[6];
				cov(1,0) = startNodeInfo[6];  
				cov(0,2) = startNodeInfo[7];
				cov(2,0) = startNodeInfo[7];
				cov(1,1) = startNodeInfo[8]; 
				cov(1,2) = startNodeInfo[9]; 
				cov(2,1) = startNodeInfo[9]; 
				cov(2,2) = startNodeInfo[10];

				mid_vector3[0] = startNodeInfo[2];
				mid_vector3[1] = startNodeInfo[3];
				mid_vector3[2] = startNodeInfo[4];

				nextTrans.fromVector(mid_vector3);	
				nextVarMatrix = cov.inverse();
					
				Jacobian_4_edge_propagate(result_synthetic_odo.first, nextTrans, j1, j2);
				covariance_propagate(result_synthetic_odo.second, nextVarMatrix, j1, j2, midVarMatrix);
				result_synthetic_odo.second = midVarMatrix;
				midTrans = (result_synthetic_odo.first) * nextTrans;
				result_synthetic_odo.first =  midTrans;

				length[0] = length[0] + sqrt(startNodeInfo[2]*startNodeInfo[2] + startNodeInfo[3]*startNodeInfo[3] );
				length[3] = result_synthetic_odo.second(0,0);
				length[4] = result_synthetic_odo.second(1,1);
		    }
		    if(vitualStart != end)
		    {
				cout<<"the final node does not equal to end vertex"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    
		    fultileBit = judgeVarianceAndDistance_Small(length);
		}

		// cout<<"end: "<<end<<endl;
	}

	//input: the start and end node serial number correspond to the segment you want to synthesize
	//output: the tranform and covariance info of the synthesized segment, in Trans and cov
	void synthesize_odo_edges(int start, int end, std::vector<std::array<double,11>> & OdoInf, g2o::SE2 & Trans, 
		Matrix3d & cov, std::array<double, 5> & length, bool & fultileBit)
	{
		fultileBit = 0;
		std::array<double,11> startNodeInfo = OdoInf[start];
		g2o::Vector3  mid_vector3;
		Matrix3d m_m, m2, J1, J2;
		g2o::SE2 edge2;
		int lengthNode = end -start;
		if(startNodeInfo[0] != start){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(lengthNode < 0){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		cov = Matrix3d::Identity(); 
		if(lengthNode == 0){
			cov(0,0)= 0; cov(1,1)= 0; cov(2,2) = 0;
			mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
			Trans.fromVector(mid_vector3);
		}
		else
		{
			if(lengthNode == 1 and (startNodeInfo[1] != end))
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}			
			cov(0,0) = startNodeInfo[5]; 
			cov(0,1) = startNodeInfo[6];
			cov(1,0) = startNodeInfo[6];  
			cov(0,2) = startNodeInfo[7];
			cov(2,0) = startNodeInfo[7];
			cov(1,1) = startNodeInfo[8]; 
			cov(1,2) = startNodeInfo[9]; 
			cov(2,1) = startNodeInfo[9]; 
			cov(2,2) = startNodeInfo[10];

			Matrix3d cov_mid;
			cov_mid = cov.inverse();
			cov = cov_mid;

			mid_vector3[0] = startNodeInfo[2];
			mid_vector3[1] = startNodeInfo[3];
			mid_vector3[2] = startNodeInfo[4];
			Trans.fromVector(mid_vector3);

			length[0] = sqrt(Trans[0]*Trans[0]+Trans[1]*Trans[1]);
			length[1] = abs(Trans[0]);
			length[2] = abs(Trans[1]);

			length[3] = cov(0,0);
			length[4] = cov(1,1);

			int j=start, minus = start;
			for(j=start+1; j<end; j++)
			{
				m2 = Matrix3d::Identity();
				m2(0,0 ) = OdoInf[j][5];
				m2(0,1 ) = OdoInf[j][6];
				m2(0,2 ) = OdoInf[j][7];
				m2(1,0 ) = m2(0,1 );
				m2(2,0 ) = m2(0,2 );

				m2(1,1) = OdoInf[j][8]; 
				m2(1,2) = OdoInf[j][9]; 
				m2(2,1 ) = m2(1,2 );
				m2(2,2) = OdoInf[j][10];

				cov_mid = m2.inverse();
				m2 = cov_mid;

				mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
				edge2.fromVector(mid_vector3);
				length[0] = length[0] + sqrt(edge2[0]*edge2[0]+edge2[1]*edge2[1]);
				length[1] = length[1]+abs(edge2[0]);
				length[2] = length[2]+abs(edge2[1]);
				Jacobian_4_edge_propagate(Trans, edge2, J1, J2);
				covariance_propagate(cov, m2, J1, J2, m_m);
				cov = m_m;

				length[3] = cov(0,0);
				length[4] = cov(1,1);

				Trans *= edge2;//update transform

				if(j - minus == 30)
				{
					fultileBit = judgeVarianceAndDistance_Small(length);
					if(fultileBit == 1)
						return;
				}
			}
			fultileBit = judgeVarianceAndDistance_Small(length);

			if(OdoInf[j-1][1] != end)
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}
	}

	void Jacobian_4_edge_propagate(g2o::SE2 & TransA, g2o::SE2 & TransB, Matrix3d & J1, Matrix3d & J2)
	{
		J1 = Matrix3d::Identity();
		J2 = Matrix3d::Identity();
		J1(0,2) = -TransB[0]*sin(TransA[2])-TransB[1]*cos(TransA[2]);
		J1(1,2) = TransB[0]*cos(TransA[2])-TransB[1]*sin(TransA[2]);
		J2(0,0) = cos(TransA[2]);J2(0,1) = -sin(TransA[2]);J2(1,0) = sin(TransA[2]);J2(1,1) = cos(TransA[2]);
	}

	void covariance_propagate(Matrix3d & cov1, Matrix3d & cov2, Matrix3d & J1, Matrix3d & J2, Matrix3d & result)
	{
		result = J1*cov1*(J1.transpose())+J2*cov2*(J2.transpose());
	}

	std::pair<bool, double> check_single_loop(int loop_to_check, cluster & _clustersFoundi, int loop_checked,
		std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > &transSequence_cluster, int clus_num, int group_num)//, double& statis
	{
		g2o::SE2 loop1, edgeGo, loop2, edgeBack, Edge_midd, Edge_go_midd, Edge_back_midd, transform_interator;
		Matrix3d Cov_loop1, Cov_edgeGo, Cov_loop2, Cov_edgeBack, Cov_midd, Cov_go_midd, Cov_back_midd, J1, J2, Cov_interator;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		// std::array<int,8> nodes;
		int go_odo_num = 0, back_odo_num = 0;
		// //nodes loop check
		// nodes[0] = _clustersFoundi.positionserial[loop_checked][0];
		// nodes[1] = _clustersFoundi.positionserial[loop_checked][3];
		if(loop_to_check <= loop_checked)
		{
			cout<<"loop_to_check: "<<loop_to_check<<" loop_checked: "<<loop_checked<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		loop1        =  transSequence_cluster[2][loop_checked].first;
		Cov_loop1    =  transSequence_cluster[2][loop_checked].second;
		loop2        =  transSequence_cluster[2][loop_to_check].first;
		Cov_loop2    =  transSequence_cluster[2][loop_to_check].second;
		if(loop_to_check == loop_checked+1)
		{ 
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	
		}
		else
		{
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	

			for(int sumodo = 1; sumodo < loop_to_check; sumodo++)
			{

				Matrix3d m2 = Matrix3d::Identity(), m_m, J1_2_odo, J2_2_odo;

				Edge_go_midd  =  transSequence_cluster[1][loop_checked+sumodo].first;
				Cov_go_midd   =  transSequence_cluster[1][loop_checked+sumodo].second;	
				Edge_back_midd  =  transSequence_cluster[0][loop_checked+sumodo].first;
				Cov_back_midd   =  transSequence_cluster[0][loop_checked+sumodo].second;	

				Jacobian_4_edge_propagate(edgeGo, Edge_go_midd, J1, J2);
				covariance_propagate(Cov_edgeGo, Cov_go_midd, J1, J2, m_m);
				Cov_edgeGo = m_m;
				edgeGo *= Edge_go_midd;//update transform

				Jacobian_4_edge_propagate(edgeBack, Edge_back_midd, J1, J2);
				covariance_propagate(Cov_edgeBack, Cov_back_midd, J1, J2, m_m);
				Cov_edgeBack = m_m;
				edgeBack *= Edge_go_midd;//update transform

				if(loop_checked+sumodo+1 == loop_to_check)
					break;
				else if(sumodo == loop_to_check-1)
				{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}
		Jacobian_4_edge_propagate(loop1, edgeGo, J1, J2);//generate jacobian 
		covariance_propagate(Cov_loop1, Cov_edgeGo, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator = loop1 * edgeGo;//update transform

		Edge_midd = loop2.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_loop2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		Edge_midd = edgeBack.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_edgeBack, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov_interator.inverse();
		T(0)= transform_interator[0];
		T(1) = transform_interator[1];
		T(2) = transform_interator[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(1) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(2) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);
		cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		cout<<"cov: "<<endl<<Cov_interator<<endl;
		cout<<"transformDistance: "<<transformDistance<<endl;
		cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		cout<<" "<<endl;
		returnV.second = transformDistance;
		if (transformDistance < 7.81)
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;
	}	

//input is four pair<transform, covariance> of two loop and two odo edge segments
//return is one pair<pass_check_or_not, transfrom_distance>
	std::pair<bool, double> check_single_loop_inter(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter, 
	double& covX, double & covY, Matrix3d & displayCov,double & transX_residual, double & transY_residual, double & transA_residual,
	double & length, bool & futileb)//, double& statis
	{
		g2o::SE2 loop1, loop2, Edge_midd;
		Matrix3d Cov1, Cov2, Cov_midd, J1, J2;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		loop1        =  transSequence_cluster_inter[0].first;
		Cov1    =  transSequence_cluster_inter[0].second;

		for(int i = 1; i<4; i++)
		{
			loop2   =  transSequence_cluster_inter[i].first;
			Cov2    =  transSequence_cluster_inter[i].second;

			Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
			covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
			Cov1 = Cov_midd;
			// //dispaly the covariacne matrix
			// cout<<Cov2.inverse()<<endl;
			// cout<<" "<<endl;

			loop1 = loop1 * loop2;//update transform
		}
		// exit(0);

		// double SNR = (length*length)/sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1));
		double SNR = (length)/sqrt(Cov1(0, 0)*Cov1(0, 0) + Cov1(1, 1)*Cov1(1, 1));

		// cout<<"length[0]: "<<length[0]<<" length[3]: "<<length[3]<<" length[4]: "<<length[4]<<endl;
		// cout<<"sqrt(length[3]*length[3] + length[4]*length[4]): "<<sqrt(length[3]*length[3] + length[4]*length[4])<<endl;
		if(SNR < snrThres)//3.16 corresponding to 5 dB
		// if(SNR < 1)
		// if(SNR < 10)
		{
			futileb = 1;
		}
		else
			futileb = 0;


		displayCov = Cov1;
		covX = Cov1(0, 0);
		covY = Cov1(1, 1);
		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov1.inverse();
		T(0)= loop1[0];
		T(1) = loop1[1];
		T(2) = loop1[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		transX_residual = T(0);
		transY_residual = T(1);
		transA_residual = T(2);
		// cout<<"transX_residual: "<<transX_residual<<"transY_residual: "<<transY_residual<<"transA_residual: "<<transA_residual<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// exit(0);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);

		// cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		// cout<<"cov: "<<endl<<Cov_interator<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		// cout<<" "<<endl;
		returnV.second = transformDistance;
		if (abs(transformDistance) < 7.815)//11.34 correspond to P=0.01, which 7.81 to 0.05
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}

	void inter_cons_check(std::vector<std::vector<int> > & consistent_pair_clusterr_real, 
		std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> >  > > &transSequence_cluster)
	{
		
	}

	IntPairSet& getClusterByID(int id){
		return clusterIDtoLoopsMap[id];
	}

	size_t clusterCount()
	{
		return _clustersFound.size();
	}

	bool deleteCluster(int clusterID)
	{
		clusterIDtoLoopsMap.erase(clusterID);

		for(IntPairIDMap::iterator it= loopToClusterIDMap.begin();
				it!=loopToClusterIDMap.end(); it++)
		{
			if(it->second == clusterID)
				loopToClusterIDMap.erase(it->first);
		}

		return true;
	}
	void merge_cluster(std::vector<std::vector<int> > & consistent_pair_clusterr_real)
	{
	}
		
};

#endif /* CLUSTER_HPP_ */
