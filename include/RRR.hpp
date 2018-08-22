// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef RRR_HPP_
#define RRR_HPP_


#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream> 
#include <sys/time.h>

#include "g2o_Interface.hpp"
#include "cluster.hpp"
#include "utils.hpp"

#include "g2o_interface_se2_zihao.hpp"


template <class BackEndInterface> //!< @param BackEndInterace one of the back-end
class RRR
{
	Clusterizer clusterizer;
	

	double clusteringThreshold;
	int nIterations;
	int ID_IGNORE;
	string cluster_filename;
	string OptimizedResultWithoutRecheck;//file name of the optimized result with     out     doubt loop check
	string OptimizedResultWithRecheck;	//file name of the optimized result with doubt loop check
	string OptimizedResultWithReRecheck, OptimizedResultBeforInter;
	string final_suvived_loops_fileName, OptimizedResultBeforactive, OptimizedResultAfteractive;
	string OptimizedResultIntra ;
	string OptimizedResultIntraInter ;
	string OptimizedResultIntraConflict ;
	string OptimizedResultIntraConflictInter ;	
   // const char *OPre=OptiRefile.data();
	std::vector<std::pair<std::pair<int, int>, double> >  odoEdgeRelateLC_Error;
	std::vector<std::array<double, 5>>  odoEdgeError;
	int edgeDimension;

	IntSet _goodSet;
	// double fiveNineBelief = 28*2;
	double twoNineBelief = 11.345, fiveNineBelief = 28, nineNineBelief = 44.8413, twelveBelief = 58.9413, upperLimitChi = 77.4, tenUpperLimit = 774;


public:
	std::vector<std::pair<int, double>> SmallSumChiCluster;
	std::pair<int, double> eleForSamllSCCluster;
	find_element find_ele;
	BackEndInterface* gWrapper;
	std::map<int, double> overallChi2SET4singleElementCluster;
	int small_numIterations = 3;

	RRR(
		double ClusteringThreshold = 1,
		int numIterations = 5,
		string Cluster_filename = "unmane"
			):
	clusteringThreshold(ClusteringThreshold),
				nIterations(numIterations),
				cluster_filename(Cluster_filename)
	{
		ID_IGNORE 			= -2;
		gWrapper 			= new BackEndInterface;
		edgeDimension 		= gWrapper->edgeDimension();

		string fullName = cluster_filename;
		int posfile = fullName.find_last_of('/');
		string path = fullName.substr(0, posfile+1); 
		string  fileName=fullName.substr(posfile+1);
	    fileName = fileName.substr(0, fileName.length()-4);

	    // resultfile = path+"result_"+fileName+'_'+".g2o";

		// clusterizer.nameofclusterfile   =    "zihao_cluster"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
		clusterizer.nameofclusterfile   =    "zihao_cluster"+fileName+"_.g2o";
		final_suvived_loops_fileName = "final_suvived_loops_"+fileName+".g2o";

		OptimizedResultWithoutRecheck = path+"result_zihao_"+fileName+"_.g2o";
		OptimizedResultWithRecheck = path+"result_zihao_"+fileName+"_recheck_.g2o";
		OptimizedResultWithReRecheck = path+"result_zihao_"+fileName+"_rerecheck_.g2o";
		OptimizedResultBeforInter =  path+"result_zihao_"+fileName+".g2o";
		OptimizedResultBeforactive=  path+"result_zihao_"+fileName+"_ba.g2o";
		OptimizedResultAfteractive=  path+"result_zihao_"+fileName+"_aa.g2o";

		OptimizedResultIntra =  path+"result_zihao_"+fileName+"_1.g2o";
		OptimizedResultIntraInter = path+"result_zihao_"+fileName+"_2.g2o";
		OptimizedResultIntraConflict = path+"result_zihao_"+fileName+"_3.g2o";
		OptimizedResultIntraConflictInter = path+"result_zihao_"+fileName+"_4.g2o";	
	}

	bool chi2_test(const std::vector<double> & dis_clu,std::vector<double> & eachChi2Error)
	{
		eachChi2Error.clear();
		double p_value ;
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		double eachError = 0.0;
		for(int i = 0; i < dis_clu.size(); i++)
		{
			eachError = (dis_clu[i]-mean)*(dis_clu[i]-mean);
			eachChi2Error.push_back(eachError);
		}
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value = 1-chi2_p(dis_clu.size()-1, chi_statis);	

		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		if(p_value > 0.05)
			return 1;
		else
			return 0; 
	}

	bool read(const char* filename)
	{
		// cout<<" read re why exit quietly? "<<endl;
		IntPairSet loops;
		bool status = gWrapper->read(filename);		// Read the g2o file
		if(!status) return status;
		gWrapper->getLoopClosures(loops);
		
		std::cerr<<"Number of loop closures : "<<loops.size()<<std::endl;
		
		//clusterizer.clusterize(loops,clusteringThreshold);	// Make clusters
		clusterizer.clusterize_zihao(loops,filename);
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		return true;
	}

	bool setOptimizer(void* optimizerPtr, const char* filename)
	{
		cout<<"xx re why exit quietly? "<<endl;
		IntPairSet loops;
		gWrapper->setOptimizer(optimizerPtr);
		gWrapper->getLoopClosures(loops);
		
		std::cerr<<"Number of loop closures : "<<loops.size()<<std::endl;

		// clusterizer.clusterize(loops,clusteringThreshold);	// Make clusters
		clusterizer.clusterize_zihao(loops,filename);
	
		return true;
	}

	bool intraClusterConsistent_second(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop)
	{
		cout<<" "<<endl;
		cout<<"second is running"<<endl;
		// std::cerr<<" chi2 %95 2: "<<utils::chi2(2)<<std::endl;
		IntPairDoubleMap chi2LinkErrors, chi2LinkErrors2;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);
		int originalTotalNum=0;
		originalTotalNum=currentCluster.size();
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<double>  eachChi2Error, dis_clu;
		std::vector<std::pair<std::pair<int,int> , double> > chi2AndLoopPair;
		std::pair<std::pair<int,int> , double> eleffssxzh;
		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
		chiStstis4eachLoop.clear();
		doubtLoopNumberVector.clear();

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		// if( currentCluster.find(std::pair<int,int>(1875, 0)) != currentCluster.end())

		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		double meanError = 0;
		int countLoop = 0;
		std::vector<std::pair<int,int> > badloop ;
		IntPairSet good_set, good_set2, good_set3;//bad_set
		std::vector<std::pair<int,int> > ID_IGNORE_vector;

		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0;
		std::pair<int, int> biggestErrorLoop;
		biggestErrorLoop.first = -1;
		biggestErrorLoop.second = -1;

		cout<<" "<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			eleffssxzh.first  = eIt->first;
			eleffssxzh.second = eIt->second;
			chi2AndLoopPair.push_back(eleffssxzh);
			dis_clu.push_back(eIt->second);

			cout<<"eIt->first.first: "<<eIt->first.first<<" and "<<eIt->first.second<<" 's cluster id is: "<<clusterizer.loopToClusterIDMap[eIt->first]<<endl;

			// if(eIt->first.first == 9016 and eIt->first.second == 3876)
			// {
			// 	cout<<"get the loop's id when add to cluster "<<clusterizer.loopToClusterIDMap[eIt->first]<<endl;
			// 	std::cin.get();
			// }

			// if(eIt->second > biggestErr)
			// // if(eIt->second > 100)
			// {
			// 	biggestErr = eIt->second;
				
			// 	biggestErrorLoop.first  = eIt->first.first;
			// 	biggestErrorLoop.second = eIt->first.second;
			// }
		}
		// if(clusterID == 85 and clusterizer.loopToClusterIDMap[std::pair<int,int> (9016, 3876)] != 85)
		// 	exit(0);

		if(biggestErr > upperLimitChi and originalTotalNum==1)
		{
			cout<<"biggestErr: "<<biggestErr<<endl;
			return 0;
		}

		//delete the loop withgiggest error that bigger than upperlimit then redo opeimize		
		if(biggestErr > upperLimitChi )
		{

			clusterizer.setClusterID(biggestErrorLoop,ID_IGNORE);
			currentCluster.clear();
			currentCluster = clusterizer.getClusterByID(clusterID);
			originalTotalNum = currentCluster.size();
			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		}

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		cout<<"*********************** first time to find the bad loops ************************"<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
			utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			
			chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
			meanError = meanError + eIt->second;
			countLoop = countLoop+1;

			if(eIt->second > utils::chi2(edgeDimension))
			// if(eIt->second > 100)
			{
				// bad_set.insert(eIt->first);
				badloop.push_back(eIt->first);//save the bad loop nodes
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

			}
			else
			{
				cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				good_set.insert(eIt->first);
			}
		}
		cout<<"**********************************************************************************"<<endl;
		cout<<" "<<endl;
		// eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

		cout<<" "<<endl;
		cout<<"*********************** optimise on good loop set 1 ************************"<<endl;
		if(good_set.size() > 0 and badloop.size() > 0 )
		{
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(good_set, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
				utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				
				chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
				meanError = meanError + eIt->second;
				countLoop = countLoop+1;

				if(eIt->second > utils::chi2(edgeDimension))
				// if(eIt->second > 100)
				{
					// bad_set.insert(eIt->first);
					badloop.push_back(eIt->first);//save the bad loop nodes
					cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

				}
				else
				{
					cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
					good_set2.insert(eIt->first);
				}
			}
			cout<<"**********************************************************************************"<<endl;
			cout<<" "<<endl;

			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
			int goodLoopNum1 = good_set.size(),goodLoopNum2 = good_set2.size();

			while(goodLoopNum1 > goodLoopNum2)
			{
				goodLoopNum1 = goodLoopNum2;
				good_set3.clear();
				cout<<" "<<endl;
				cout<<"*********************** optimise on good loop set 2 ************************"<<endl;

					chi2LinkErrors.clear();
					chi2LinkErrors = gWrapper->optimize(good_set2, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
					activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
					activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

				eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

				cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
					utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
				for(; eIt!=eEnd; eIt++)
				{
					if(eIt->first.first < 0)
						continue;
					
					chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
					meanError = meanError + eIt->second;
					countLoop = countLoop+1;

					if(eIt->second > utils::chi2(edgeDimension))
					// if(eIt->second > 100)
					{
						// bad_set.insert(eIt->first);
						badloop.push_back(eIt->first);//save the bad loop nodes
						cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

					}
					else
					{
						cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
						good_set3.insert(eIt->first);
					}
				}
				cout<<"**********************************************************************************"<<endl;
				cout<<" "<<endl;
				goodLoopNum2 = good_set3.size();
				if(goodLoopNum2 == 0)
					break;
				if(goodLoopNum2 == goodLoopNum1)
					break;
				good_set2.clear();
				good_set2.insert(good_set3.begin(), good_set3.end());
			}
		}



		//although badloop equal to original elements number ,we shoudl not delete it without futher check
		if(badloop.size() == originalTotalNum)
		{
			cout<<"the bad loop number equals to the total loop number so delete the cluster"<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}
		// //this is an bad if condition
		// if(1 == originalTotalNum and activeChi2Graph > utils::chi2_continuous(edgeDimension* ceil(activeEdgeCount/10), 0.1))
		// {
		// 	cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
		// 	cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<endl;

		// 	cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
		// 	std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;

		// 	return false;
		// } 
		
		if((originalTotalNum > 1) and (badloop.size() > 0) and (originalTotalNum - badloop.size())> 0)
		{
			
			//16.266 correspond to 0.001
			//6.251                0.1
			if(badloop.size() == 1)
			{
				if(chi2LinkErrors[badloop[0]] < upperLimitChi)
					doubtLoopNumberVector.push_back(badloop[0]);
				else
					ID_IGNORE_vector.push_back(badloop[0]);
				//else judge it as bad one
			}
			else
			{
				IntPairDoubleMap re_chi2LinkErrors;
				for(int i = 0; i < badloop.size(); i++)
				{
					IntPairSet mix_set_to_reintracheck;
					std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
					//put loops in good set in to mixset 
					for(; recheck != endrecheck; recheck++){
						mix_set_to_reintracheck.insert(*recheck);
					}
					//put a single bad loop into the mixset 
					mix_set_to_reintracheck.insert(badloop[i]);
					re_chi2LinkErrors.clear();
					re_chi2LinkErrors = gWrapper->optimize(mix_set_to_reintracheck, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
					cout<<"loop "<<badloop[i].first<<" "<<badloop[i].second <<" has a residual error as "<<re_chi2LinkErrors[badloop[i]]<<endl; 
					// cout<<"chi 9999 : "<<chi2_continuous(3, 0.9999)<<endl; 
					// cout<<"chi 99999 : "<<chi2_continuous(3, 0.99999)<<endl; 
					// cout<<"chi 999999 : "<<chi2_continuous(3, 0.999999)<<endl; 
					if(re_chi2LinkErrors[badloop[i]] > upperLimitChi or 
						re_chi2LinkErrors[IntPair(-2,0)] > utils::chi2_continuous(re_chi2LinkErrors[IntPair(-2,-1)]* 3, 0.95)
						or re_chi2LinkErrors[IntPair(-1,0)] > utils::chi2_continuous(re_chi2LinkErrors[IntPair(-1,-1)]* 3, 0.95))//100//16.266 
					{
						//it is a bad one 
						ID_IGNORE_vector.push_back(badloop[i]);	
					}
					else if(re_chi2LinkErrors[badloop[i]] > utils::chi2(edgeDimension))
					{
						doubtLoopNumberVector.push_back(badloop[i]);
					}
					else
						good_set.insert(badloop[i]);
					//
				}
			}
			//redo intracheck, first clear the cluster, then add good set ele into it, which is followed by adding doubt loop into it	
			currentCluster.clear();
			std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
			//add goodset element into it
			for(; recheck != endrecheck; recheck++)
			{
				currentCluster.insert(*recheck);
			}
			//add doubt element into it
			for(int i = 0; i < doubtLoopNumberVector.size(); i++)
			{
				currentCluster.insert(doubtLoopNumberVector[i]);
			}
			//redo optimize
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
		}

		cout<<"currentCluster.size: "<<currentCluster.size()<<endl;

		activeChi2Graph = chi2LinkErrors[IntPair(-1, 0)];
		activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		//clear the doubtset
		doubtLoopNumberVector.clear();

		//if the overall residual if less then the limit
	
		if( activeChi2Graph < utils::chi2(edgeDimension* activeEdgeCount) and 
			chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95))
		{
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			//pick out the loops whose residual error is larger than the limit
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;

				// if(eIt->second > utils::chi2(edgeDimension))
				if(eIt->second > upperLimitChi)
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is deleted"<<endl;
					clusterizer.setClusterID(eIt->first,ID_IGNORE);
					chiStatisVector.push_back(eIt->second);
				}
				else if(eIt->second > utils::chi2(edgeDimension))
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is put into double vector"<<endl;

					chiStatisVector.push_back(eIt->second);
					doubtLoopNumberVector.push_back(eIt->first); 
					clusterizer.setClusterID(eIt->first,ID_IGNORE);
				}
			}
			//put the loops failed in first check into ID_IGNORE
			for(int i = 0; i < ID_IGNORE_vector.size(); i++)
			{
				clusterizer.setClusterID(ID_IGNORE_vector[i],ID_IGNORE);
			}
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;

			if((clusterizer.getClusterByID(clusterID).size() + doubtLoopNumberVector.size())/double(realOrigianlNum) < 0.49)
			{
				cout<<"bad loops account for a propertion that is larger than 0.5 , so we delete this cluster"<<endl;
				sleep(5);
				return 0;
			}
			if(clusterizer.getClusterByID(clusterID).size() >= 1 )
			{
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				if(originalTotalNum == 1)
				{
					overallChi2SET4singleElementCluster[clusterID] = activeChi2Graph;	
				}
				return true;
			}
			else
			{
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
				return false;
			}
		}
		//if the overall error is bigger than the limit, the original RRR delete the cluster without other opertation,but it needs some change
		//first, if the number of loops is bigger than one pick out the smallest residual error to judge whether less than limit
		//second, 
		else
		{
			//if the total chi statis is beyond the limit, then find
			cout<<" "<<endl;
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			
			cout<<"sum chi2 error of loops: "<<chi2LinkErrors[IntPair(-2,0)]<<" loop number: "<<chi2LinkErrors[IntPair(-2,-1)]<<
				" "<<utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95)<<endl;

			std::cerr<<" Cluster "<<clusterID<<" has been eliminated!"<<std::endl;
			return false;
		}
		// exit(0);

		return false;
	}
bool intraClusterConsistent_second_active(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop, string ba, string aa)
	{
		cout<<" "<<endl;
		cout<<"active is running"<<endl;

		
		std::vector<double>  eachChi2Error, dis_clu;
		std::vector<std::pair<std::pair<int,int> , double> > chi2AndLoopPair;
		std::pair<std::pair<int,int> , double> eleffssxzh;

		// std::cerr<<" chi2 %95 2: "<<utils::chi2(2)<<std::endl;
		IntPairDoubleMap chi2LinkErrors, chi2LinkErrors2;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);

		int originalTotalNum = currentCluster.size();
		originalLoopNum = originalTotalNum;

		int realOrigianlNum = originalLoopNum;

		//get the node range
		std::set<int> clusterNodeSet;
		int startRange = 0, endRange = 0;

						// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);

		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
		for(; findRange != findRangeEnd; findRange++)
		{
			clusterNodeSet.insert((*findRange).first);
			clusterNodeSet.insert((*findRange).second);
		}
		startRange = *clusterNodeSet.begin();
		endRange   = *clusterNodeSet.rbegin();

		//pick out the subgraph corresponding to the node range
		chi2LinkErrors =  gWrapper->optimize_active(currentCluster, small_numIterations, startRange, endRange, odoEdgeError, ba, aa);
		// chi2LinkErrors =  gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

						// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
		// loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();
		// loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		// loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();

		chiStstis4eachLoop.clear();
		doubtLoopNumberVector.clear();

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		// int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		int   activeEdgeCount = endRange- startRange + originalTotalNum;

		// if( currentCluster.find(std::pair<int,int>(1875, 0)) != currentCluster.end())
		double meanError = 0;
		int countLoop = 0;
		std::vector<std::pair<int,int> > badloop ;
		IntPairSet good_set, good_set2, good_set3;//bad_set
		std::vector<std::pair<int,int> > ID_IGNORE_vector;

		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr, smallestErr;
		std::pair<int, int> biggestErrorLoop;
		biggestErrorLoop.first = -1;
		biggestErrorLoop.second = -1;

		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		cout<<"*********************** first time to find the bad loops ************************"<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
			utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			eleffssxzh.first  = eIt->first;
			eleffssxzh.second = eIt->second;
			chi2AndLoopPair.push_back(eleffssxzh);
			dis_clu.push_back(eIt->second);

			chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
			meanError = meanError + eIt->second;
			countLoop = countLoop+1;

			if(eIt->second > utils::chi2(edgeDimension))
			// if(eIt->second > 100)
			{
				// bad_set.insert(eIt->first);
				badloop.push_back(eIt->first);//save the bad loop nodes
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
			}
			else
			{
				cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				good_set.insert(eIt->first);
			}
			cout<<"cluster id is: "<<clusterizer.loopToClusterIDMap[eIt->first]<<endl;
		}
		// if(clusterID == 85 and clusterizer.loopToClusterIDMap[std::pair<int,int> (9016, 3876)] != 85)
		// 	exit(0);
		biggestErr  = *(std::max_element(std::begin(dis_clu), std::end(dis_clu)));
		smallestErr = *(std::min_element(std::begin(dis_clu), std::end(dis_clu)));

		// if(chi2LinkErrors[IntPair(-1,0)] > utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-1,-1)], 0.95))  
		//activeEdgeCount = endRange- startRange + originalTotalNum;
		if(chi2LinkErrors[IntPair(-1,0)] > utils::chi2_continuous(edgeDimension * activeEdgeCount, 0.95))
		{
			cout<<"whole graph chi2 error: "<<chi2LinkErrors[IntPair(-1,0)] <<" is bigger than "<<utils::chi2_continuous(edgeDimension * activeEdgeCount, 0.95)<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return 0;
		}

		if(smallestErr > upperLimitChi)
		{
			cout<<"abandon this cluster because smallestErr: "<<smallestErr<<" is bigger than upperlimit"<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return 0;
		}

		if(chi2LinkErrors[IntPair(-2,0)] > utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95))
		{
			cout<<"chi2 error on all loop closures: "<<smallestErr<<" is bigger than "<<utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95)<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return 0;
		}
		cout<<"**********************************************************************************"<<endl;
		// eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		if(good_set.size() > 0 and badloop.size() > 0 )
		{
			cout<<" "<<endl;
			cout<<"*********************** optimise on good loop set ************************"<<endl;

			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(good_set, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = good_set.size()+ endRange - startRange;// chi2LinkErrors[IntPair(-1,-1)];

			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
				utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
				meanError = meanError + eIt->second;
				countLoop = countLoop+1;

				if(eIt->second > utils::chi2(edgeDimension))
				// if(eIt->second > 100)
				{
					// bad_set.insert(eIt->first);
					badloop.push_back(eIt->first);//save the bad loop nodes
					cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				}
				else
				{
					cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
					good_set2.insert(eIt->first);
				}
			}
			cout<<"**********************************************************************************"<<endl;
			cout<<" "<<endl;

			int goodLoopNum1 = good_set.size(),goodLoopNum2 = good_set2.size();

			if(goodLoopNum1 != goodLoopNum2)
			{
				cout<<"goodSET has bad loops, so we delete it"<<endl;
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
				return 0;
			}
		}

		//although badloop equal to original elements number ,we shoudl not delete it without futher check
		if(badloop.size() == originalTotalNum)
		{
			cout<<"the bad loop number equals to the total loop number so delete the cluster"<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}
		// //this is an bad if condition
		// if(1 == originalTotalNum and activeChi2Graph > utils::chi2_continuous(edgeDimension* ceil(activeEdgeCount/10), 0.1))
		// {
		// 	cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
		// 	cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<endl;

		// 	cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
		// 	std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;

		// 	return false;
		// } 
		
		if((originalTotalNum > 1) and (badloop.size() > 1) and (originalTotalNum - badloop.size())> 0)
		{
			
			//16.266 correspond to 0.001
			//6.251                0.1
			if(badloop.size() == 1)
			{
				if(chi2LinkErrors[badloop[0]] < upperLimitChi)
					doubtLoopNumberVector.push_back(badloop[0]);
				else
					ID_IGNORE_vector.push_back(badloop[0]);
				//else judge it as bad one
			}
			else
			{
				IntPairDoubleMap re_chi2LinkErrors;
				IntPairSet mix_set_to_reintracheck;
				for(int i = 0; i < badloop.size(); i++)
				{
					mix_set_to_reintracheck.clear();
					std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
					//put loops in good set in to mixset 
					for(; recheck != endrecheck; recheck++){
						mix_set_to_reintracheck.insert(*recheck);
					}
					//put a single bad loop into the mixset 
					mix_set_to_reintracheck.insert(badloop[i]);
					re_chi2LinkErrors.clear();
					re_chi2LinkErrors = gWrapper->optimize(mix_set_to_reintracheck, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
					cout<<"loop "<<badloop[i].first<<" "<<badloop[i].second <<" has a residual error as "<<re_chi2LinkErrors[badloop[i]]<<endl; 
					// cout<<"chi 9999 : "<<chi2_continuous(3, 0.9999)<<endl; 
					// cout<<"chi 99999 : "<<chi2_continuous(3, 0.99999)<<endl; 
					// cout<<"chi 999999 : "<<chi2_continuous(3, 0.999999)<<endl; 
					if(re_chi2LinkErrors[badloop[i]] > upperLimitChi or 
						re_chi2LinkErrors[IntPair(-2,0)] > utils::chi2_continuous((mix_set_to_reintracheck.size()+ endRange - startRange)* 3, 0.95)
						or re_chi2LinkErrors[IntPair(-1,0)] > utils::chi2_continuous(re_chi2LinkErrors[IntPair(-1,-1)]* 3, 0.95))//100//16.266 
					{
						//it is a bad one 
						ID_IGNORE_vector.push_back(badloop[i]);	
					}
					else if(re_chi2LinkErrors[badloop[i]] > utils::chi2(edgeDimension))
					{
						doubtLoopNumberVector.push_back(badloop[i]);
					}
					else
						good_set.insert(badloop[i]);
					//
				}
			}
			//redo intracheck, first clear the cluster, then add good set ele into it, which is followed by adding doubt loop into it	
			currentCluster.clear();
			std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
			//add goodset element into it
			for(; recheck != endrecheck; recheck++)
			{
				currentCluster.insert(*recheck);
			}
			//add doubt element into it
			for(int i = 0; i < doubtLoopNumberVector.size(); i++)
			{
				currentCluster.insert(doubtLoopNumberVector[i]);
			}
			//redo optimize
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
		}

		cout<<"currentCluster.size: "<<currentCluster.size()<<endl;

		activeChi2Graph = chi2LinkErrors[IntPair(-1, 0)];
		// activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		//clear the doubtset
		doubtLoopNumberVector.clear();

		//if the overall residual if less then the limit
	
		if( activeChi2Graph < utils::chi2(edgeDimension* activeEdgeCount) and 
			chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95))
		{
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			//pick out the loops whose residual error is larger than the limit
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;

				// if(eIt->second > utils::chi2(edgeDimension))
				if(eIt->second > upperLimitChi)
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is deleted"<<endl;
					// clusterizer.setClusterID(eIt->first,ID_IGNORE);
					chiStatisVector.push_back(eIt->second);
				}
				else if(eIt->second > utils::chi2(edgeDimension))
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is put into double vector"<<endl;

					chiStatisVector.push_back(eIt->second);
					doubtLoopNumberVector.push_back(eIt->first); 
					// clusterizer.setClusterID(eIt->first,ID_IGNORE);
				}
			}
			//put the loops failed in first check into ID_IGNORE
			for(int i = 0; i < ID_IGNORE_vector.size(); i++)
			{
				// clusterizer.setClusterID(ID_IGNORE_vector[i],ID_IGNORE);
			}
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;

			if((clusterizer.getClusterByID(clusterID).size() + doubtLoopNumberVector.size())/double(realOrigianlNum) < 0.49)
			{
				cout<<"bad loops account for a propertion that is larger than 0.5 , so we delete this cluster"<<endl;
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
				sleep(5);
				return 0;
			}
			if(clusterizer.getClusterByID(clusterID).size() >= 1 )
			{
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				if(originalTotalNum == 1)
				{
					overallChi2SET4singleElementCluster[clusterID] = activeChi2Graph;	
				}
				return true;
			}
			else
			{
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
				return false;
			}
		}
		//if the overall error is bigger than the limit, the original RRR delete the cluster without other opertation,but it needs some change
		//first, if the number of loops is bigger than one pick out the smallest residual error to judge whether less than limit
		//second, 
		else
		{
			//if the total chi statis is beyond the limit, then find
			cout<<" "<<endl;
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			
			cout<<"sum chi2 error of loops: "<<chi2LinkErrors[IntPair(-2,0)]<<" loop number: "<<chi2LinkErrors[IntPair(-2,-1)]<<
				" "<<utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95)<<endl;

			std::cerr<<" Cluster "<<clusterID<<" has been eliminated!"<<std::endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}
		// exit(0);

		return false;
	}
	//chiStstis4eachLoop  is vector to store the chis and nodes of loops. 
	bool intraClusterConsistent(int clusterID, std::vector<double> & chiStatisVector, std::vector<std::pair<int, int> > & doubtLoopNumberVector, 
		int & originalLoopNum, std::vector<std::pair<double, std::pair<int, int> > > & chiStstis4eachLoop)
	{
		// std::cerr<<" chi2 %95 2: "<<utils::chi2(2)<<std::endl;
		IntPairDoubleMap chi2LinkErrors, chi2LinkErrors2;
		IntPairSet& currentCluster = clusterizer.getClusterByID(clusterID);
		int originalTotalNum=0;
		originalTotalNum=currentCluster.size();
		originalLoopNum = originalTotalNum;
		int realOrigianlNum = originalLoopNum;
		std::vector<double>  eachChi2Error, dis_clu;
		std::vector<std::pair<std::pair<int,int> , double> > chi2AndLoopPair;
		std::pair<std::pair<int,int> , double> eleffssxzh;
		//the third parameter should be set to negative things, it is used to skip some certain loop when debug
		chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
		chiStstis4eachLoop.clear();
		doubtLoopNumberVector.clear();

		float activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		int   activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		// if( currentCluster.find(std::pair<int,int>(1875, 0)) != currentCluster.end())

		if(originalTotalNum == 1)
		{		
			for(int i = 0; i < odoEdgeError.size(); i++)
			{
				// allerr1=allerr1+ odoEdgeError[i][0];
				// allerr2=allerr2+ odoEdgeError1[i][0];
				// angleError = odoEdgeError[i][3]*odoEdgeError[i][3]*odoEdgeError[i][4];
				// cout<<odoEdgeError[i][0]<<" "<<odoEdgeError[i][4] <<endl;//" "<<odoEdgeError[i][3]*odoEdgeError[i][3]*1420618.101850<<" "
					// << angleError<<endl;

				// if(odoEdgeError[i][4]> 0.11)//0.25
				// {
				// 	cout<<"residual on single edge: "<<odoEdgeError[i][4]<<endl;
				// 	cout<<"it is one member cluster "<<clusterID<<endl;
				// 	sleep(5);
				// 	return false;
				// }
			}
		}

		IntPairDoubleMap::iterator eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		double meanError = 0;
		int countLoop = 0;
		std::vector<std::pair<int,int> > badloop ;
		IntPairSet good_set;//bad_set
		std::vector<std::pair<int,int> > ID_IGNORE_vector;

		//check if there is someone loop whose chi2 error is bigger than the signgle edge limit error
		double biggestErr = 0;
		std::pair<int, int> biggestErrorLoop;
		biggestErrorLoop.first = -1;
		biggestErrorLoop.second = -1;

		cout<<" "<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			eleffssxzh.first  = eIt->first;
			eleffssxzh.second = eIt->second;
			chi2AndLoopPair.push_back(eleffssxzh);
			dis_clu.push_back(eIt->second);
			// if(eIt->second > biggestErr)
			// // if(eIt->second > 100)
			// {
			// 	biggestErr = eIt->second;
				
			// 	biggestErrorLoop.first  = eIt->first.first;
			// 	biggestErrorLoop.second = eIt->first.second;
			// }
		}

		if(biggestErr > upperLimitChi and originalTotalNum==1)
		{
			cout<<"biggestErr: "<<biggestErr<<endl;
			return 0;
		}

		//delete the loop withgiggest error that bigger than upperlimit then redo opeimize		
		if(biggestErr > upperLimitChi )
		{

			clusterizer.setClusterID(biggestErrorLoop,ID_IGNORE);
			currentCluster.clear();
			currentCluster = clusterizer.getClusterByID(clusterID);
			originalTotalNum = currentCluster.size();
			//the third parameter should be set to negative things, it is used to skip some certain loop when debug
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
			activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];
		}

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		cout<<" "<<endl;
		cout<<"*********************** first time to find the bad loops ************************"<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<
			utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
		for(; eIt!=eEnd; eIt++)
		{
			if(eIt->first.first < 0)
				continue;
			
			chiStstis4eachLoop.push_back(std::pair<double, std::pair<int, int> > (eIt->second, eIt->first));
			meanError = meanError + eIt->second;
			countLoop = countLoop+1;
	
			if(eIt->second > utils::chi2(edgeDimension))
			// if(eIt->second > 100)
			{
				// bad_set.insert(eIt->first);
				badloop.push_back(eIt->first);//save the bad loop nodes
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;

			}
			else
			{
				cout<<"good loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				good_set.insert(eIt->first);
			}
		}
		cout<<"**********************************************************************************"<<endl;
		cout<<" "<<endl;



		//although badloop equal to original elements number ,we shoudl not delete it without futher check
		if(badloop.size() == originalTotalNum)
		{
			cout<<"the bad loop number equals to the total loop number so delete the cluster"<<endl;
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
			return false;
		}
		//this is an bad if condition
		if(1 == originalTotalNum and activeChi2Graph > utils::chi2_continuous(edgeDimension* ceil(activeEdgeCount/10), 0.1))
		{
			cout<<"this cluster has one loop, and the activeChi2Graph is bigger than the overall limit, so delete this cluster"<<endl;
			cout<<"activeChi2Graph: "<<activeChi2Graph <<" activeEdgeCount: "<<activeEdgeCount<<endl;

			cout<<"gWrapper->odoSize: "<<gWrapper->odoSize<<endl;		
			std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;

			return false;
		}

		// //sort chiStatisWhenAdd2Cluster
		// std::sort(clusterizer._clustersFound[clusterID].chiStatisWhenAdd2Cluster.begin(), 
		// clusterizer._clustersFound[clusterID].chiStatisWhenAdd2Cluster.end(), cmp);  //排序 

		// std::sort(chiStstis4eachLoop.begin(), chiStstis4eachLoop.end(), cmp);  //排序 
		
		if((originalTotalNum > 1) and (badloop.size() > 0) and (originalTotalNum - badloop.size())> 0)
		{
			
			//16.266 correspond to 0.001
			//6.251                0.1
			if(badloop.size() == 1)
			{
				if(chi2LinkErrors[badloop[0]] < upperLimitChi)
					doubtLoopNumberVector.push_back(badloop[0]);
				else
					ID_IGNORE_vector.push_back(badloop[0]);
				//else judge it as bad one
			}
			else
			{
				IntPairDoubleMap re_chi2LinkErrors;
				for(int i = 0; i < badloop.size(); i++)
				{
					IntPairSet mix_set_to_reintracheck;
					std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
					//put loops in good set in to mixset 
					for(; recheck != endrecheck; recheck++){
						mix_set_to_reintracheck.insert(*recheck);
					}
					//put a single bad loop into the mixset 
					mix_set_to_reintracheck.insert(badloop[i]);
					re_chi2LinkErrors.clear();
					re_chi2LinkErrors = gWrapper->optimize(mix_set_to_reintracheck, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
					cout<<"previous bad loop "<<badloop[i].first<<" "<<badloop[i].second <<" has a residual error as "<<re_chi2LinkErrors[badloop[i]]<<endl; 
				
					recheck= good_set.begin();//endrecheck = good_set.end();
					//put loops in good set in to mixset 
					for(; recheck != endrecheck; recheck++){
						cout<<"good set's loop  "<<(*recheck).first<<" "<<(*recheck).second <<" has a residual error as "<<re_chi2LinkErrors[*recheck]<<endl; 
					}
					// cout<<"chi 9999 : "<<chi2_continuous(3, 0.9999)<<endl; 
					// cout<<"chi 99999 : "<<chi2_continuous(3, 0.99999)<<endl; 
					// cout<<"chi 999999 : "<<chi2_continuous(3, 0.999999)<<endl; 
					if(re_chi2LinkErrors[badloop[i]] > upperLimitChi or 
						re_chi2LinkErrors[IntPair(-2,0)] > utils::chi2_continuous(re_chi2LinkErrors[IntPair(-2,-1)]* 3, 0.95)
						or re_chi2LinkErrors[IntPair(-1,0)] > utils::chi2_continuous(re_chi2LinkErrors[IntPair(-1,-1)]* 3, 0.95))//100//16.266 
					{
						//it is a bad one 
						ID_IGNORE_vector.push_back(badloop[i]);	
					}
					else if(re_chi2LinkErrors[badloop[i]] > utils::chi2(edgeDimension))
					{
						doubtLoopNumberVector.push_back(badloop[i]);
					}
					else
						good_set.insert(badloop[i]);
					//
				}

			}
			// std::cin.get();

			//redo intracheck, first clear the cluster, then add good set ele into it, which is followed by adding doubt loop into it	
			currentCluster.clear();
			std::set<std::pair<int,int> >::iterator recheck= good_set.begin(), endrecheck = good_set.end();
			//add goodset element into it
			for(; recheck != endrecheck; recheck++)
			{
				currentCluster.insert(*recheck);
			}
			//add doubt element into it
			for(int i = 0; i < doubtLoopNumberVector.size(); i++)
			{
				currentCluster.insert(doubtLoopNumberVector[i]);
			}
			//redo optimize
			chi2LinkErrors.clear();
			chi2LinkErrors = gWrapper->optimize(currentCluster, small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
		}

		cout<<"currentCluster.size: "<<currentCluster.size()<<endl;

		activeChi2Graph = chi2LinkErrors[IntPair(-1,0)];
		activeEdgeCount = chi2LinkErrors[IntPair(-1,-1)];

		eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();
		//clear the doubtset
		doubtLoopNumberVector.clear();

		//if the overall residual if less then the limit
	
		if( activeChi2Graph < utils::chi2(edgeDimension* activeEdgeCount) and 
			chi2LinkErrors[IntPair(-2,0)] < utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95))
		{
			eIt = chi2LinkErrors.begin(), eEnd=chi2LinkErrors.end();

			//pick out the loops whose residual error is larger than the limit
			for(; eIt!=eEnd; eIt++)
			{
				if(eIt->first.first < 0)
					continue;
				cout<<"loop "<<eIt->first.first<<" "<<eIt->first.second<<" has residual error "<<eIt->second<<endl;
				// if(eIt->second > utils::chi2(edgeDimension))
				if(eIt->second > upperLimitChi)
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is deleted"<<endl;
					clusterizer.setClusterID(eIt->first,ID_IGNORE);
					chiStatisVector.push_back(eIt->second);
				}
				else if(eIt->second > utils::chi2(edgeDimension))
				{
					//std::cerr<<" edgeDimension: "<<edgeDimension<<std::endl;
					cout<<eIt->first.first<<" "<<eIt->first.second<<" is put into double vector"<<endl;

					chiStatisVector.push_back(eIt->second);
					doubtLoopNumberVector.push_back(eIt->first); 
					clusterizer.setClusterID(eIt->first,ID_IGNORE);
				}
			}

			cout<<"one intracheck end"<<endl;
			// std::cin.get();

			//put the loops failed in first check into ID_IGNORE
			for(int i = 0; i < ID_IGNORE_vector.size(); i++)
			{
				clusterizer.setClusterID(ID_IGNORE_vector[i],ID_IGNORE);
			}
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;

			if((clusterizer.getClusterByID(clusterID).size() + doubtLoopNumberVector.size())/double(realOrigianlNum) < 0.49)
			{
				cout<<"bad loops account for a propertion that is larger than 0.5 , so we delete this cluster"<<endl;
				sleep(5);
				return 0;
			}
			if(clusterizer.getClusterByID(clusterID).size() >= 1 )
			{
				std::cerr<<" Cluster "<<clusterID <<"survived with "<<clusterizer.getClusterByID(clusterID).size()<<" links "<<"from "<<realOrigianlNum<<std::endl;
				if(originalTotalNum == 1)
				{
					overallChi2SET4singleElementCluster[clusterID] = activeChi2Graph;	
				}
				return true;
			}
			else
			{
				std::cerr<<" Empty Cluster "<<clusterID<<" has  "<<originalTotalNum<<" links"<<std::endl;
				return false;
			}
		}
		//if the overall error is bigger than the limit, the original RRR delete the cluster without other opertation,but it needs some change
		//first, if the number of loops is bigger than one pick out the smallest residual error to judge whether less than limit
		//second, 
		else
		{
			//if the total chi statis is beyond the limit, then find
			cout<<" "<<endl;
			cout<<"activeChi2Graph: "<<activeChi2Graph << " utils::chi2: "<<utils::chi2(edgeDimension* activeEdgeCount)<<" activeEdgeCount: "<<activeEdgeCount<<endl;
			
			cout<<"sum chi2 error of loops: "<<chi2LinkErrors[IntPair(-2,0)]<<" loop number: "<<chi2LinkErrors[IntPair(-2,-1)]<<
				" "<<utils::chi2_continuous(edgeDimension * chi2LinkErrors[IntPair(-2,-1)], 0.95)<<endl;

			std::cerr<<" Cluster "<<clusterID<<" has been eliminated!"<<std::endl;
			return false;
		}
		// exit(0);

		return false;
	}

	
	// eigen: edge_propagate_se2(g2o::VertexSE2* e1, g2o::VertexSE2* e2)
	// {
	// 	g2o::Vector3 p = (dynamic_cast<g2o::VertexSE2*>(rrr.gWrapper->optimizer->vertex(1000)))->estimate().toVector();
	// }


	bool gatherLinks(const IntSet& clusterList, IntPairSet& links)
	{
		IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
		// cout<<"gather links form cluster: ";
		for( ; it!=end; it++)
		{
			// cout<<(*it)<<" ";
			IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
			links.insert(currentCluster.begin(),currentCluster.end());
		}
		// cout<<endl;
		return true;
	}
	bool gatherLinks(const std::vector<int> & clusterList, IntPairSet& links)
	{
		// cout<<"gather links form cluster: ";
		for(int gather = 0 ; gather < clusterList.size(); gather++)
		{
			// cout<<(clusterList[gather])<<" ";
			IntPairSet& currentCluster = clusterizer.getClusterByID(clusterList[gather]);
			links.insert(currentCluster.begin(),currentCluster.end());
		}
		// cout<<endl;
		return true;
	}
	bool gatherLinks_equalElement(const IntSet& clusterList, IntPairSet& links)
	{
		IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
		int minElementNum = clusterizer.getClusterByID(*it).size();
		it++;

		for( ; it!=end; it++)
		{
			if(clusterizer.getClusterByID(*it).size() < minElementNum)
				minElementNum= clusterizer.getClusterByID(*it).size();
		}

		it = clusterList.begin(), end = clusterList.end();
		for( ; it!=end; it++)
		{
			IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
			auto pointer = currentCluster.begin();
			for(int toINsert = 0; toINsert < minElementNum; toINsert++)
			{
				if(pointer == currentCluster.end())
				{
					printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				links.insert(*(pointer));
				pointer++;
			}
		}
		return true;
	}

	double getNumberIn(const IntSet& clusterList)
	{
		IntSet::const_iterator it = clusterList.begin(), end = clusterList.end();
		double NumOfGoodLC=0;

		for( ; it!=end; it++)
		{
			IntPairSet& currentCluster = clusterizer.getClusterByID(*it);
			NumOfGoodLC += currentCluster.size();
		}
		return NumOfGoodLC;
	}	

	bool interClusterConsistent(IntSet& H, IntSet& goodSet, IntSet& RejectSet)
	{

		if(H.empty()) return true;

		IntPairSet activeLinks;

		gatherLinks(goodSet,activeLinks);
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];

		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}



		if(
				activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size())
			)
		{
			goodSet.insert(H.begin(),H.end());
			return true; // all done .. we found a consistent solution
		}
		else
		{
			// Find which cluster is causing the problems
			// Iterate over everything
			// Find the error for each cluster
			// Sum cluster-wise

			IntPairSet::iterator it = activeLinks.begin(), end = activeLinks.end();

			std::map< int, double > errorMap;

			for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
					it!=end;
					it++)
			{
				if(it->first.first >=0)
				{
					int thisLinkClusterID = clusterizer.getClusterID(it->first);
					errorMap[thisLinkClusterID]	+= it->second;
				}
			}

			IntSet::iterator sIt = H.begin(), sEnd = H.end();

			double min_CI = 0;
			int rejectID = -1;

			for ( ; sIt!=sEnd; sIt++)
			{
				double CI = errorMap[*sIt]/utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size());

				if( CI >= min_CI ) // Just looking for the ones in H
				{
					min_CI = CI;
					rejectID = *sIt;
				}
			}
			H.erase(rejectID);
			RejectSet.insert(rejectID);

			return interClusterConsistent(H,goodSet,RejectSet);
		}
		return true;
	}

	bool interClusterConsistent_second(IntSet& H, IntSet& goodSet, IntSet& RejectSet)
	{
		if(H.empty()) return true;

		IntPairSet activeLinks;

		gatherLinks(goodSet,activeLinks);
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];

		double allLinksError = linkErrors[IntPair(-2,0)];
		if(activeLinks.size() != linkErrors[IntPair(-2,-1)])
		{
				cout<<"activeLinks.size(): "<<activeLinks.size()<<" linkErrors[IntPair(-2,0)]: "<<linkErrors[IntPair(-2,0)]<<endl;
				printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		}
		// for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
		// 		it!=end;
		// 		it++)
		// {
		// 	if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
		// 		allLinksError += it->second;
		// }



		if(
				activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size())
			)
		{
			goodSet.insert(H.begin(),H.end());
			return true; // all done .. we found a consistent solution
		}
		else
		{
			// Find which cluster is causing the problems
			// Iterate over everything
			// Find the error for each cluster
			// Sum cluster-wise

			IntPairSet::iterator it = activeLinks.begin(), end = activeLinks.end();

			std::map< int, double > errorMap;

			for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
					it!=end;
					it++)
			{
				if(it->first.first >=0)
				{
					int thisLinkClusterID = clusterizer.getClusterID(it->first);
					errorMap[thisLinkClusterID]	+= it->second;
				}
			}

			IntSet::iterator sIt = H.begin(), sEnd = H.end();

			double min_CI = 0;
			int rejectID = -1;

			for ( ; sIt!=sEnd; sIt++)
			{
				double CI = errorMap[*sIt]/utils::chi2(edgeDimension*clusterizer.getClusterByID(*sIt).size());

				if( CI >= min_CI ) // Just looking for the ones in H
				{
					min_CI = CI;
					rejectID = *sIt;
				}
			}
			cout<<"delete cluster "<<rejectID<<endl;
			cout<<"its mean chi2 error is "<<min_CI<<endl;
			cout<<"activeChi2Graph: "<<activeChi2Graph<< "  threshold: "<<utils::chi2(edgeDimension*activeEdgeCount)<<endl;
			cout<<"allLinksError:   "<<allLinksError<< "  threshold: "<<utils::chi2(edgeDimension*activeLinks.size())<<endl;

			H.erase(rejectID);
			RejectSet.insert(rejectID);

			return interClusterConsistent_second(H,goodSet,RejectSet);
		}
		return true;
	}

	
bool multiClusterCheck(std::vector<int>  & A)
	{
		std::set<int> H;
		H.insert(A.begin(), A.end());
		IntPairSet activeLinks;
		std::vector<std::pair<std::pair<int,int>, double> > num;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		// if(H.size() != 2)
		// {
		// 	cout<<"two cluster check H.size: "<<H.size()<<endl;
		// 	cout<<"doubt check H should has two elements"<<endl;
		// 	exit(0);
		// }
		std::pair<std::pair<int, int> , double> ele;
		cout<<"clusters: "<<*it;
		for(; it != eit; it++)
		{


			ele.first.first = *it;
			ele.first.second = 0;
			ele.second = 0.0;
			num.push_back(ele);
		}
		cout<<endl;
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{
				// cout<<"the loop corresponding to node "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
				// 	clusterizer.getClusterID(it->first)<<" has a residual error as: "<<it->second<<endl;
				int find = 0;
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
				for(int ff1617 = 0; ff1617 < num.size(); ff1617++)
				{
					if(num[ff1617].first.first == thisLinkClusterID)
					{
						num[ff1617].first.second++;
						num[ff1617].second += it->second;
						find = 1;
						break;
					}
				}
				if(find == 0)
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"can not find the cluster to add the chi2 error"<<endl;
					exit(0);
				}
			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}

		int everyclusterPass = 1;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			cout<<"cluster: "<<num[clusterNum].first.first<<" has "<<num[clusterNum].first.second<<" loops"<<endl;
			cout<<"its chi2 error is "<<num[clusterNum].second<<" with chi2 threshold "<<
				utils::chi2(edgeDimension*(num[clusterNum].first.second))<<endl;
			if(num[clusterNum].second >  utils::chi2(edgeDimension*(num[clusterNum].first.second)))
			{
				cout<<"fail to set everyclusterPass as 1"<<endl;
				everyclusterPass = 0;
			}
		}
		if(activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size()
				and everyclusterPass))
		{
			return true; // all done .. we found a consistent solution
		}
		else
			return false;
	}

bool multiClusterCheck(std::set<int> H)
	{

		IntPairSet activeLinks;
		std::vector<std::pair<std::pair<int,int>, double> > num;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		// if(H.size() != 2)
		// {
		// 	cout<<"two cluster check H.size: "<<H.size()<<endl;
		// 	cout<<"doubt check H should has two elements"<<endl;
		// 	exit(0);
		// }
		std::pair<std::pair<int, int> , double> ele;
		cout<<"clusters: "<<*it;
		for(; it != eit; it++)
		{


			ele.first.first = *it;
			ele.first.second = 0;
			ele.second = 0.0;
			num.push_back(ele);
		}
		cout<<endl;
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{
				// cout<<"the loop corresponding to node "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
				// 	clusterizer.getClusterID(it->first)<<" has a residual error as: "<<it->second<<endl;
				int find = 0;
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
				for(int ff1617 = 0; ff1617 < num.size(); ff1617++)
				{
					if(num[ff1617].first.first == thisLinkClusterID)
					{
						num[ff1617].first.second++;
						num[ff1617].second += it->second;
						find = 1;
						break;
					}
				}
				if(find == 0)
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}

		int everyclusterPass = 1;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			cout<<"cluster: "<<num[clusterNum].first.first<<" has "<<num[clusterNum].first.second<<" loops"<<endl;
			cout<<"its chi2 error is "<<num[clusterNum].second<<" with chi2 threshold "<<
				utils::chi2(edgeDimension*(num[clusterNum].first.second))<<endl;
			if(num[clusterNum].second >  utils::chi2(edgeDimension*(num[clusterNum].first.second)))
			{
				cout<<"fail to set everyclusterPass as 1"<<endl;
				everyclusterPass = 0;
			}
		}
		if(activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size()
				and everyclusterPass))
		{
			return true; // all done .. we found a consistent solution
		}
		else
			return false;
	}

bool multiClusterCheck(std::set<int> & H, std::set<int> & H_refine)
	{
		H_refine.clear();
		IntPairSet activeLinks;
		std::vector<std::pair<std::pair<int,int>, double> > num;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		// if(H.size() != 2)
		// {
		// 	cout<<"two cluster check H.size: "<<H.size()<<endl;
		// 	cout<<"doubt check H should has two elements"<<endl;
		// 	exit(0);
		// }
		cout<<" "<<endl;
		std::pair<std::pair<int, int> , double> ele;
		cout<<"clusters: ";
		for(; it != eit; it++)
		{

			cout<<*it<<" ";
			ele.first.first = *it;
			ele.first.second = 0;
			ele.second = 0.0;
			num.push_back(ele);
		}
		cout<<endl;
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{

				// cout<<"the loop corresponding to node "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
				// 	clusterizer.getClusterID(it->first)<<" has a residual error as: "<<it->second<<endl;
				int find = 0;
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
				for(int ff1617 = 0; ff1617 < num.size(); ff1617++)
				{
					if(num[ff1617].first.first == thisLinkClusterID)
					{
						cout<<"chi2error is "<<it->second<<" in cluster "<<thisLinkClusterID<<endl;
						num[ff1617].first.second++;
						num[ff1617].second += it->second;
						find = 1;
						break;
					}
				}

				if(find == 0)
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}
		cout<<"allLinksError  : "<<allLinksError<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph<<endl;
		int everyclusterPass = 1;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			cout<<"cluster: "<<num[clusterNum].first.first<<" has "<<num[clusterNum].first.second<<" loops"<<endl;
			cout<<"its chi2 error is "<<num[clusterNum].second<<" with chi2 threshold "<<utils::chi2(edgeDimension*(num[clusterNum].first.second))<<endl;
			if(num[clusterNum].second >  utils::chi2(edgeDimension*(num[clusterNum].first.second)))
			{
				cout<<"fail to set everyclusterPass as 1"<<endl;
				everyclusterPass = 0;
				// break;
			}
			else
				H_refine.insert(num[clusterNum].first.first);
			// everyclusterPass = 1;
		}

		if(activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size()
				and
				everyclusterPass)
			)
		{
			return true; // all done .. we found a consistent solution
		}
		else
		{
			return false;
		}
	}

bool multiClusterCheck(std::vector<int> A, std::set<int> & H_refine)
	{
		std::set<int> H;
		H_refine.clear();
		H.insert(A.begin(), A.end());
		IntPairSet activeLinks;
		std::vector<std::pair<std::pair<int,int>, double> > num;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		// if(H.size() != 2)
		// {
		// 	cout<<"two cluster check H.size: "<<H.size()<<endl;
		// 	cout<<"doubt check H should has two elements"<<endl;
		// 	exit(0);
		// }
		cout<<" "<<endl;
		std::pair<std::pair<int, int> , double> ele;
		cout<<"clusters: ";
		for(; it != eit; it++)
		{

			cout<<*it<<" ";
			ele.first.first = *it;
			ele.first.second = 0;
			ele.second = 0.0;
			num.push_back(ele);
		}
		cout<<endl;
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{

				// cout<<"the loop corresponding to node "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
				// 	clusterizer.getClusterID(it->first)<<" has a residual error as: "<<it->second<<endl;
				int find = 0;
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
				for(int ff1617 = 0; ff1617 < num.size(); ff1617++)
				{
					if(num[ff1617].first.first == thisLinkClusterID)
					{
						cout<<"chi2error is "<<it->second<<" in cluster "<<thisLinkClusterID<<endl;
						num[ff1617].first.second++;
						num[ff1617].second += it->second;
						find = 1;
						break;
					}
				}

				if(find == 0)
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;

					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}
		cout<<"allLinksError  : "<<allLinksError<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph<<endl;
		int everyclusterPass = 1;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			cout<<"cluster: "<<num[clusterNum].first.first<<" has "<<num[clusterNum].first.second<<" loops"<<endl;
			cout<<"its chi2 error is "<<num[clusterNum].second<<" with chi2 threshold "<<utils::chi2(edgeDimension*(num[clusterNum].first.second))<<endl;
			if(num[clusterNum].second >  utils::chi2(edgeDimension*(num[clusterNum].first.second)))
			{
				cout<<"fail to set everyclusterPass as 1"<<endl;
				everyclusterPass = 0;
				// break;
			}
			else
				H_refine.insert(num[clusterNum].first.first);
			// everyclusterPass = 1;
		}

		if(activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size()
				and
				everyclusterPass)
			)
		{
			return true; // all done .. we found a consistent solution
		}
		else
		{
			return false;
		}
	}

	bool multiClusterCheck_onlyMinusOnePerTime(std::set<int> & H, std::set<int> & H_refine)
	{
		// std::set<int> H;
		// H_refine.clear();
		// H.insert(A.begin(), A.end());
		IntPairSet activeLinks;
		std::vector<std::pair<std::pair<int,int>, double> > num;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		// if(H.size() != 2)
		// {
		// 	cout<<"two cluster check H.size: "<<H.size()<<endl;
		// 	cout<<"doubt check H should has two elements"<<endl;
		// 	exit(0);
		// }
		cout<<" "<<endl;
		std::pair<std::pair<int, int> , double> ele;
		cout<<"clusters: ";
		for(; it != eit; it++)
		{

			cout<<*it<<" ";
			ele.first.first = *it;
			ele.first.second = 0;
			ele.second = 0.0;
			num.push_back(ele);
		}
		cout<<endl;
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{

				// cout<<"the loop corresponding to node "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
				// 	clusterizer.getClusterID(it->first)<<" has a residual error as: "<<it->second<<endl;
				int find = 0;
				int thisLinkClusterID = clusterizer.getClusterID(it->first);
				errorMap[thisLinkClusterID]	+= it->second;
				for(int ff1617 = 0; ff1617 < num.size(); ff1617++)
				{
					if(num[ff1617].first.first == thisLinkClusterID)
					{
						cout<<"chi2error is "<<it->second<<" in cluster "<<thisLinkClusterID<<endl;
						num[ff1617].first.second++;
						num[ff1617].second += it->second;
						find = 1;
						break;
					}
				}

				if(find == 0)
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}
		cout<<"allLinksError  : "<<allLinksError<<endl;
		cout<<"activeChi2Graph: "<<activeChi2Graph<<endl;

		std::vector<std::pair<int, double>> infoLargeChi2ErrCluster;
		std::pair<int, double> lele;
		int everyclusterPass = 1;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			cout<<"cluster: "<<num[clusterNum].first.first<<" has "<<num[clusterNum].first.second<<" loops"<<endl;
			cout<<"its chi2 error is "<<num[clusterNum].second<<" with chi2 threshold "<<utils::chi2(edgeDimension*(num[clusterNum].first.second))<<endl;
			if(num[clusterNum].second >  utils::chi2(edgeDimension*(num[clusterNum].first.second)))
			{
				cout<<"fail to set everyclusterPass as 1"<<endl;
				lele.first  = num[clusterNum].first.first;
				lele.second = (num[clusterNum].second)/num[clusterNum].first.second;
				infoLargeChi2ErrCluster.push_back(lele);				
				everyclusterPass = 0;
				// break;
			}
			else
			{
				// H_refine.insert(num[clusterNum].first.first);
			}
			// everyclusterPass = 1;
		}

		int sessseei= 0;
		double dsfnw =0.0;
		if(everyclusterPass == 0)
		{
			for(int sf9 = 0; sf9 < infoLargeChi2ErrCluster.size(); sf9++)
			{
				if(infoLargeChi2ErrCluster[sf9].second > dsfnw)
				{
					dsfnw = infoLargeChi2ErrCluster[sf9].second;
					sessseei = infoLargeChi2ErrCluster[sf9].first;
				}
			}
		}
		cout<<"delete cluster "<<sessseei<<endl;
		for(int clusterNum = 0; clusterNum < num.size(); clusterNum++)
		{
			if(num[clusterNum].first.first == sessseei)
				continue;
			else
			{
				H_refine.insert(num[clusterNum].first.first);
			}
			// everyclusterPass = 1;
		}


		if(activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size()
				and
				everyclusterPass)
			)
		{
			return true; // all done .. we found a consistent solution
		}
		else
		{
			return false;
		}
	}

bool twoClusterCheck(std::set<int> H)
	{
		int thisLinkClusterID;
		IntPairSet activeLinks;
		std::pair<std::pair<int,double>, int> num1,num2;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		if(H.size() != 2)
		{
			cout<<"two cluster check H.size: "<<H.size()<<endl;
			cout<<"doubt check H should has two elements"<<endl;
			exit(0);
		}
		int wuna = 0;
		for(; it != eit; it++)
		{
			cout<<"H element: "<<*it<<endl;
			if(wuna == 0)
				num1.first.first = *it;
			else
				num2.first.first = *it;
			wuna = wuna+1;
		}
		num1.first.second = 0;
		num2.first.second = 0;	
		num1.second = 0;
		num2.second = 0;	
		gatherLinks(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{

				cout<<"loop "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
					clusterizer.getClusterID(it->first)<<" has an error as: "<<it->second<<endl;
		
				thisLinkClusterID = clusterizer.getClusterID(it->first);

				errorMap[thisLinkClusterID]	+= it->second;
				if(num1.first.first == thisLinkClusterID)
				{
					num1.first.second = num1.first.second + it->second;
					num1.second = num1.second + 1;

				}
				else if(num2.first.first == thisLinkClusterID)
				{
					num2.first.second = num2.first.second + it->second;
					num2.second = num2.second + 1;
				}
				else
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"it->first: "<<it->first.first<<" "<<it->first.second<<endl;
					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}			
		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}

		if(
				activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size())
				and
				num2.first.second   <   utils::chi2(edgeDimension*(num2.second))
				and
				num1.first.second   <   utils::chi2(edgeDimension*(num1.second))
			)
		{
			return true; // all done .. we found a consistent solution
		}
		else
		{
			cout<<"activeChi2Graph: "<<activeChi2Graph<<" < "<<utils::chi2(edgeDimension*activeEdgeCount)<<endl;
			cout<<"allLinksError: "<<allLinksError<<" < "<<utils::chi2(edgeDimension*activeLinks.size())<<endl;
			cout<<"num2.first.second: "<<num2.first.second<<" < "<<utils::chi2(edgeDimension*(num2.second))<<endl;
			cout<<"num1.first.second: "<<num1.first.second<<" < "<<utils::chi2(edgeDimension*(num1.second))<<endl;									
			return false;
		}
	}

	bool twoClusterCheck_equalElement(std::set<int> H, int & toDelete)
	{
		int thisLinkClusterID;
		IntPairSet activeLinks;
		std::pair<std::pair<int,double>, int> num1,num2;
		std::set<int>::iterator it = H.begin(), eit = H.end();
		if(H.size() != 2)
		{
			cout<<"two cluster check H.size: "<<H.size()<<endl;
			cout<<"doubt check H should has two elements"<<endl;
			exit(0);
		}
		int wuna = 0;
		for(; it != eit; it++)
		{
			cout<<"H element: "<<*it<<endl;
			if(wuna == 0)
				num1.first.first = *it;
			else
				num2.first.first = *it;
			wuna = wuna+1;
		}
		num1.first.second = 0;
		num2.first.second = 0;	
		num1.second = 0;
		num2.second = 0;	
		gatherLinks_equalElement(H,activeLinks);

		IntPairDoubleMap
			linkErrors = gWrapper->optimize(activeLinks,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);

		std::map< int, double > errorMap;
		// cout<<" "<<endl;
		// cout<<"two cluster check between "<<num1.first.first<<" and "<<num2.first.first<<endl;
		for( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >=0)
			{

				cout<<"loop "<<it->first.first<<" "<< it->first.second<<" in cluster "<< 
					clusterizer.getClusterID(it->first)<<" has an error as: "<<it->second<<endl;
		
				thisLinkClusterID = clusterizer.getClusterID(it->first);

				errorMap[thisLinkClusterID]	+= it->second;
				if(num1.first.first == thisLinkClusterID)
				{
					num1.first.second = num1.first.second + it->second;
					num1.second = num1.second + 1;

				}
				else if(num2.first.first == thisLinkClusterID)
				{
					num2.first.second = num2.first.second + it->second;
					num2.second = num2.second + 1;
				}
				else
				{
					cout<<"thisLinkClusterID : "<<thisLinkClusterID<<endl;
					cout<<"it->first: "<<it->first.first<<" "<<it->first.second<<endl;
					cout<<"two cluster check but find the third cluster"<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}

		if(num2.second != num1.second)		
		{
			cout<<"this error is in equal TWO CLUSTER check , but two cluster number do not equal"<<endl;
			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		else
		{
			if(num2.first.second/num2.second > num1.first.second/num1.second)
			{
				// if(num2.first.second > utils::chi2(edgeDimension*(num2.second)))
				if(num2.first.second > utils::chi2_continuous(edgeDimension*(num2.second), 0.95))
					toDelete = num2.first.first;
				else
					toDelete = -3;

			}
			else
			{
				// if(num1.first.second > utils::chi2(edgeDimension*(num1.second)))
				if(num1.first.second >  utils::chi2_continuous(edgeDimension*(num2.second), 0.95))
				{
					toDelete = num1.first.first;
				}
				else
					toDelete = -3;
			}

			cout<<"the cluster to delete is "<<toDelete<<endl;
		}	


		float activeChi2Graph = linkErrors[IntPair(-1,0)];
		int   activeEdgeCount = linkErrors[IntPair(-1,-1)];
		double allLinksError = 0;
		for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
				it!=end;
				it++)
		{
			if(it->first.first >= 0) // Special case for the two values we returned with a -1 as the first
				allLinksError += it->second;
		}

		if(
				activeChi2Graph 	< 	utils::chi2(edgeDimension*activeEdgeCount)
				and
				allLinksError 		< 	utils::chi2(edgeDimension*activeLinks.size())
				and
				num2.first.second   <   utils::chi2(edgeDimension*(num2.second))
				and
				num1.first.second   <   utils::chi2(edgeDimension*(num1.second))
			)
		{
			return true; // all done .. we found a consistent solution
		}
		else
		{
			cout<<"activeChi2Graph: "<<activeChi2Graph<<" < "<<utils::chi2(edgeDimension*activeEdgeCount)<<endl;
			cout<<"allLinksError: "<<allLinksError<<" < "<<utils::chi2(edgeDimension*activeLinks.size())<<endl;
			cout<<"num2.first.second: "<<num2.first.second<<" < "<<utils::chi2(edgeDimension*(num2.second))<<endl;
			cout<<"num1.first.second: "<<num1.first.second<<" < "<<utils::chi2(edgeDimension*(num1.second))<<endl;									
			return false;
		}
	}

	bool robustify(bool eraseIncorrectLinks=false)
	{

		if(!gWrapper->isInitialized())
		{
			std::cerr<<" Please read in a graph file with read() or set the optimizer with setOptimizer()"<<std::endl;
			return false;
		}
		// First look for intra-cluster consistency

		IntSet
			// consistentClusters,
			hypotheses,
			hypotheses_largeCluster,
			goodSet,
			rejectSet,
			tempRejectSet;
		std::vector<int> consistentClusters;
			
		std::vector<std::pair<int, int> > loopNUmbersInConsistentClusters;
		std::vector<std::pair<int,int> > reAccept;
		
		std::cout<<"Number of Clusters found : "<<clusterizer.clusterCount()<<std::endl;
		std::cout<<"Checking Intra cluster consistency : "<<std::endl;
		//save all loops information once the cluster survive
		ofstream fileStreamr; 
		fileStreamr.open("loops_info_of_all_suvived_clusters_after_intracheck.txt",ios::trunc);
		std::pair<g2o::SE2, Matrix3d> tSave;
		std::pair<int,int> ty;
		int  originalLoopNum;
		bool passIntrachenck_second, passIntrachenck_active = 1;
		std::vector<std::pair<double, std::pair<int, int> > >  chiStstis4eachLoop;
		std::vector<std::pair<int, int> >  doubtLoopVector;
		std::vector<int> clusterIDofDoubtLoops, ToDeleteFromConflict;

		std::vector<int> acitveEffectCount, secondEffectCount;


			struct timeval t1, t2;
			gettimeofday(&t1, NULL);

		for(size_t i=0; i< clusterizer.clusterCount(); i++)
		{
			std::cout<<i<<" "; std::cout.flush();
			std::vector<double>  chiStatisVector;
			std::vector<std::pair<int, int> >  doubtLoopNumberVector;

			passIntrachenck_second = intraClusterConsistent_second(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop);		
			
			if(find_ele.find(clusterizer.vectorNewSplitCloser, i) != -1)
			{
				cout<<"find an cluster that need to be split but not"<<endl;
				std::cin.get();
			}


	
			// passIntrachenck_active = intraClusterConsistent_second_active(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop,
			// 			OptimizedResultBeforactive,	OptimizedResultAfteractive);




			// currentCluster.clear();
			// currentCluster = clusterizer.getClusterByID(clusterID);

			// 		IntPairSet::iterator findRange = currentCluster.begin(), findRangeEnd = currentCluster.end();
			// for(; findRange != findRangeEnd; findRange++)
			// {
			// 	clusterNodeSet.insert((*findRange).first);
			// 	clusterNodeSet.insert((*findRange).second);
			// }
			// startRange = *clusterNodeSet.begin();
			// endRange   = *clusterNodeSet.rbegin();
			// int freeDom = endRange - startRange;




			if(passIntrachenck_active == true and passIntrachenck_second == true)
			{
				cout<<"both pass"<<endl;			
			}

			if(passIntrachenck_active != true and passIntrachenck_second != true)
			{
				cout<<"both fail"<<endl;			
			}

			if(passIntrachenck_active == true and passIntrachenck_second != true)
			{
				cout<<"fail in Intrachenck_second, pass active"<<endl;	
				secondEffectCount.push_back(i);
			}

			if(passIntrachenck_active != true and passIntrachenck_second == true)
			{
				cout<<"fail in active, pass passIntrachenck_second"<<endl;
				acitveEffectCount.push_back(i);
			}		

			// std::cin.get();

			//save the suvived clusters that pass intra check
			if(passIntrachenck_second)
			{
			// 	for(int doubt = 0; doubt < doubtLoopNumberVector.size(); doubt++)
			// 	{

			// 		doubtLoopVector.push_back(doubtLoopNumberVector[doubt]);
			// 		clusterIDofDoubtLoops.push_back(i);
			// 	}
			// 	fileStreamr<<i<<"\n";
			// 	//save loops information in i to a txt file
			// 	IntPairSet::iterator it = clusterizer.getClusterByID(i).begin(), itend = clusterizer.getClusterByID(i).end();
			// 	for(; it != itend; it++)
			// 	{
			// 		int first = (*it).first, toNode = (*it).second;
			// 		int j;
			// 		for(j = 0; j<clusterizer._clustersFound[i].positionserial.size();j++)
			// 		{
			// 			if( first == clusterizer._clustersFound[i].positionserial[j][0] and
			// 				toNode == clusterizer._clustersFound[i].positionserial[j][3] )
			// 			{
			// 				ty.first = clusterizer._clustersFound[i].positionserial[j][0];
			// 				ty.second = clusterizer._clustersFound[i].positionserial[j][3];

			// 				tSave = clusterizer.LP_Trans_Covar_Map[ty];
			// 				Matrix3d fsave;
			// 				fsave = tSave.second;

			// 				fileStreamr<<"EDGE_SE2 "<<ty.first<<" "<<ty.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
			// 					<<" "<<fsave.inverse()(0,0)<<" "<<fsave.inverse()(0,1)<<" "<<fsave.inverse()(0,2)<<" "<<
			// 					fsave.inverse()(1,1)<<" "<<fsave.inverse()(1,2)<<" "<<fsave.inverse()(2,2)<<" "<<"\n";	
			// 					break;
			// 			}
			// 			if (j == clusterizer._clustersFound[i].positionserial.size()-1)
			// 			{
			// 				cout<<"can not find the postion info from the cluster for the suvived loop "<<first<<" "<<toNode<<endl;
			// 				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			// 				exit(0);
			// 			}
	
			// 		}
			// 	}	

				// for(int cbb = 0; cbb < clusterizer._clustersFound[i].positionserial.size(); cbb++)
				// {
				// 	ty.first = clusterizer._clustersFound[i].positionserial[cbb][0];
				// 	ty.second = clusterizer._clustersFound[i].positionserial[cbb][3];

				// 	tSave = clusterizer.LP_Trans_Covar_Map[ty];
				// 	Matrix3d fsave;
				// 	fsave = tSave.second;

				// 	fileStreamr<<"EDGE_SE2 "<<ty.first<<" "<<ty.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
				// 		<<" "<<fsave.inverse()(0,0)<<" "<<fsave.inverse()(0,1)<<" "<<fsave.inverse()(0,2)<<" "<<
				// 		fsave.inverse()(1,1)<<" "<<fsave.inverse()(1,2)<<" "<<fsave.inverse()(2,2)<<" "<<"\n";
				// }
				// consistentClusters.insert(i);
				consistentClusters.push_back(i);
				loopNUmbersInConsistentClusters.push_back(std::pair<int, int> (i, originalLoopNum));
			}
		}
				std::cout<<"consistentClusters size: "<<consistentClusters.size()<<std::endl;
				// sleep(2);


			gettimeofday(&t2, NULL);
			//那么函数f运行所花的时间为
			double deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec;// 微秒
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"intra time consumed: "<<deltaT<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;

		fileStreamr.close();

		// //print the clusters that suvived
		// cout<<"suvived clusters after intra check:"<<endl;
		// // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// for(std::vector<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// {
		// 	//get the cluster serial survived
		// 	//use  all_conflict_cluster to check 
		// 	cout<<(*conflict_check_iter)<<" ";
		// }
		// cout<<endl;

		// cout<<"suvived clusters but fail to pass active:"<<endl;
		// // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// for(int to_print = 0; to_print < acitveEffectCount.size(); to_print++)
		// {
		// 	//use  all_conflict_cluster to check 
		// 	cout<<acitveEffectCount[to_print]<<" ";
		// }
		// if(acitveEffectCount.size() == 0)
		// 	cout<<"no cluster fail to pass active check";
		// cout<<endl;

		// cout<<"suvived clusters from active but fail to pass second:"<<endl;
		// // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// for(int to_print = 0; to_print < secondEffectCount.size(); to_print++)
		// {
		// 	//use  all_conflict_cluster to check 
		// 	cout<<secondEffectCount[to_print]<<" ";
		// }
		// if(secondEffectCount.size() == 0)
		// 	cout<<"no cluster pass active but fail to pass second check";
		// cout<<endl;

			
		// std::cin.get();

		// //print out loops numbers  in each suvived cluster
		// cout<<"the loop numbers in suvived clusters:"<<endl;
		// for(int i = 0; i < loopNUmbersInConsistentClusters.size(); i++ )
		// { 
		// 	cout<<loopNUmbersInConsistentClusters[i].first<<"("<<loopNUmbersInConsistentClusters[i].second<<")"<<" ";
		// }
		// cout<<endl;
		// cout<<"the ACTUAL loop numbers in suvived clusters (the actual number of loops may be less because of intra check):"<<endl;
		// for(int i = 0; i < loopNUmbersInConsistentClusters.size(); i++)
		// { 
		// 	cout<<loopNUmbersInConsistentClusters[i].first<<": "<< clusterizer.getClusterByID(loopNUmbersInConsistentClusters[i].first).size()<<" "<<endl;
		// }
		// cout<<endl;
		// std::cout<<"done"<<std::endl<<std::endl;
		/// ------------- consistentCluster are self-consistent --------- ////

		//pick out the cluster whose has more then one members
		// std::set<std::set<int> > mergeSet;
		IntPairSet activeLoops;
		// IntSet::iterator
		std::vector<int>::iterator
			cIt 	= consistentClusters.begin(),
			cEnd 	= consistentClusters.end();

		// cout<<"big consistentClusters (whos loop number is bigger than 1) : "<<endl;

		for( ; cIt!=cEnd; cIt++)
		{
			// if(goodSet.find(*cIt)==goodSet.end() and
			// 		rejectSet.find(*cIt)==rejectSet.end()) // We can't ignore this because it is nether in goodSet nor in rejectSet at the moment
			// {
			// 	hypotheses.insert(*cIt);

				// if(clusterizer.getClusterByID(*cIt).size() > 1)
					hypotheses_largeCluster.insert(*cIt);
			// }
				// cout<<*cIt<<" ";
		}
		// cout<<endl;

		// std::cout<<"intra check done"<<std::endl;
		// // //pause at here 
		// std::cout << "Press \'Return\' to end." << std::endl;
		// std::cin.get();

		//print out the undertain and conflict clsuter pair
		// std::set<std::pair<int,int> > suvivedConflictClusterPairSet;


		// std::set<std::pair<int,int> > checkedConflictClusterPairSet;
		// std::pair<int,int> pairToCheck;
		std::set<int>::iterator
			cItSet	= hypotheses_largeCluster.begin(),
			cEndSet = hypotheses_largeCluster.end();
		// for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		// {
		// 	//*cIt is the serial number of the test cluster
		// 	cout<<" "<<endl;
		// 	int serialLocal = find_ele.find(clusterizer.all_uncertain_cluster, *cItSet);
		// 	cout<<"base cluster: "<< *cItSet<<endl;
		// 	// if(serialLocal != -1)
		// 	// {
		// 	// 	cout<<"uncertain clusters:"<<endl;
		// 	// 	for(int ITuncSET = 0; ITuncSET < clusterizer.all_uncertain_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
		// 	// 	{	
		// 	// 		cout<<clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]<<" ";
		// 	// 	}
		// 	// 	cout<<endl;				
		// 	// }
		// 	// else
		// 	// 	cout<<"has no uncertain cluster"<<endl;

		// 	serialLocal = find_ele.find(clusterizer.all_conflict_cluster, *cItSet);
		// 	if(serialLocal != -1)
		// 	{
				
		// 		cout<<"conflict clusters:"<<endl;
		// 		for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
		// 		{	
		// 			cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" ";

		// 		}
		// 		cout<<endl;				
		// 	}
		// 	else
		// 		cout<<"has no conflict cluster"<<endl;

		// 	if(serialLocal != -1)
		// 	{			
		// 		cout<<"suvived conflict clusters:"<<endl;
		// 		int haveSuvivedConflict;
		// 		std::pair<int,int> toFindClusterPair;
		// 		for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
		// 		{	
		// 			haveSuvivedConflict = find_ele.find(consistentClusters, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
		// 			if(haveSuvivedConflict != -1)
		// 			{
		// 				cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" "<<endl;;

		// 				// suvivedConflictClusterPairSet.insert(std::pair<int,int> ( *cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]));
		// 				// suvivedConflictClusterPairSet.insert(std::pair<int,int> ( clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET], *cItSet));

		// 				hypotheses.clear();
		// 				hypotheses.insert(*cItSet);
		// 				hypotheses.insert(clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);

		// 				if(*cItSet > clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET])
		// 				{
		// 					pairToCheck.first  = *cItSet;
		// 					pairToCheck.second = clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET];
		// 				}


		// 				//it is a map that contains conflict clsuter pair and  its distance when add to cluter
		// 				toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
		// 				if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
		// 					clusterizer.conflictClusterPairDistanceMap.end())
		// 				{
		// 					cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
		// 				}
		// 				else
		// 				{
		// 					cout<<"can not find conflict info in previous pair map"<<endl;
		// 					cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
		// 					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);

		// 					for(int ITuncSET_second = ITuncSET+1; ITuncSET_second < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET_second++)
		// 					{
		// 						toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]);
		// 						cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]<<" conflict cluster"<<endl;;

		// 						if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
		// 							clusterizer.conflictClusterPairDistanceMap.end())
		// 						{
		// 							cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
		// 						}
		// 						else
		// 						{
		// 							cout<<"can not find conflict info in previous pair map"<<endl;
		// 							cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
		// 							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 						}	
		// 					}

		// 					exit(0);
		// 				}

		// 				if(checkedConflictClusterPairSet.find(pairToCheck) != checkedConflictClusterPairSet.end())
		// 				{
		// 					cout<<"this pair has been checked already"<<endl;
		// 					continue;
		// 				}

		// 				checkedConflictClusterPairSet.insert(pairToCheck);

		// 				int  toDelete;
		// 				// if(twoClusterCheck(hypotheses))
		// 				auto passinter = twoClusterCheck_equalElement(hypotheses, toDelete);
		// 				if(passinter)
		// 				{
		// 					// ToDeleteFromConflict.push_back(toDelete);
		// 					cout<<"two confilict cluster but all pass chi2 test in optimization process"<<endl;
		// 					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 					// std::cin.get();
		// 					// exit(0);
		// 				}
		// 				else
		// 				{
		// 					if(toDelete != -3)
		// 					{
		// 						ToDeleteFromConflict.push_back(toDelete);
		// 						cout<<"the delete cluster id is "<<toDelete<<endl;
		// 						cout<<"chenk the conflict info"<<endl;
		// 						// std::cin.get();
		// 					}
		// 					else
		// 					{
		// 						cout<<"what ever"<<endl;
		// 						std::cin.get();
		// 					}
		// 				}
		// 				cout<<"cluster "<<toDelete<<" should be deleted"<<endl;
		// 			}
		// 		}
		// 		cout<<endl;				
		// 	}
		// 		//check if the uncertian cluster also exit in the survived cluster
		// }

		// cout<<"clusters to delete from conflict check are:"<<endl;
		// for(int toPrint = 0; toPrint < ToDeleteFromConflict.size(); toPrint++)
		// {
		// 	cout<<ToDeleteFromConflict[toPrint]<<" ";
		// }
		// cout<<endl;

		// cout<<"clusters suvived from intra:"<<endl;
		// // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// for(std::vector<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// {
		// 	//get the cluster serial survived
		// 	int cluster_to_test = *conflict_check_iter;
		// 	//use  all_conflict_cluster to check 
		// 	cout<<cluster_to_test<<" ";
		// }
		// cout<<endl;

		// cout<<"clusters suvived from intra and pair conflict check:"<<endl;
		// // for(std::set<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// std::vector<int> suvivedToInter;
		// for(std::vector<int>::iterator  conflict_check_iter = consistentClusters.begin(); conflict_check_iter != consistentClusters.end(); conflict_check_iter++)
		// {
		// 	//get the cluster serial survived
		// 	int cluster_to_test = *conflict_check_iter;
		// 	//use  all_conflict_cluster to check 
		// 	if(find_ele.find(ToDeleteFromConflict, cluster_to_test) == -1)
		// 	{
		// 		suvivedToInter.push_back(cluster_to_test);
		// 		cout<<cluster_to_test<<" ";
		// 	}
		// }
		// cout<<endl;

		// cout<<"all conflict suvived clsuter pairs got in clustering have been checked"<<endl;
		// cout<<" "<<endl;

		// std::cin.get();

		// std::set<int> suvivedTointerSet(suvivedToInter.begin(), suvivedToInter.end());
		// std::set<int> suvivedTointerSet_backup(suvivedToInter.begin(), suvivedToInter.end());

		std::set<int> suvivedTointerSet(consistentClusters.begin(), consistentClusters.end());
		std::set<int> suvivedTointerSet_backup(consistentClusters.begin(), consistentClusters.end());
		// std::set<int> suvivedTointerSet_backup_sec(suvivedToInter.begin(), suvivedToInter.end());

		gettimeofday(&t1, NULL);

		interClusterConsistent_second( suvivedTointerSet, goodSet, rejectSet);

					gettimeofday(&t2, NULL);
			//那么函数f运行所花的时间为
			 deltaT = (t2.tv_sec-t1.tv_sec) * 1000000 + t2.tv_usec-t1.tv_usec;// 微秒
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"inter time consumed: "<<deltaT<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;
			cout<<"*********************"<<endl;

		// interClusterConsistent_second( consistentClusters, goodSet, rejectSet);
		cout<<"number of clusters put into backend: "<<suvivedTointerSet_backup.size()<<endl;
		cout<<"number of clusters in good set: "<<goodSet.size()<<endl;


		// OptimizedResultIntra =  path+"result_zihao_"+fileName+"_Intra.g2o";
		// OptimizedResultIntraInter = path+"result_zihao_"+fileName+"_IntraInter.g2o";
		// OptimizedResultIntraConflict = path+"result_zihao_"+fileName+"_IntraConflict.g2o";
		// OptimizedResultIntraConflictInter = path+"result_zihao_"+fileName+"_IntraConflictInter.g2o";
		// string OptimizedResultIntra =  path+"result_zihao_"+fileName+"_1.g2o";
		// string OptimizedResultIntraInter = path+"result_zihao_"+fileName+"_2.g2o";
		// string OptimizedResultIntraConflict = path+"result_zihao_"+fileName+"_3.g2o";
		// string OptimizedResultIntraConflictInter = path+"result_zihao_"+fileName+"_4.g2o";		

		// cout<<"cluster set for intra and inter: ";
		// 	cItSet	= goodSet.begin();
		// 	cEndSet = goodSet.end();
		// for( ; cItSet!=cEndSet; cItSet++)
		// {
		// 	cout<<*cItSet<<" ";
		// }
		// cout<<endl;

		
		// 	cItSet	= goodSet.begin(),
		// 	cEndSet = goodSet.end();
		// for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		// {
		// 	//*cIt is the serial number of the test cluster
		// 	cout<<" "<<endl;
		// 	int serialLocal = find_ele.find(clusterizer.all_uncertain_cluster, *cItSet);
		// 	cout<<"base cluster: "<< *cItSet<<endl;

		// 	serialLocal = find_ele.find(clusterizer.all_conflict_cluster, *cItSet);
		// 	if(serialLocal != -1)
		// 	{
				
		// 		cout<<"conflict clusters:"<<endl;
		// 		for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
		// 		{	
		// 			cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" ";

		// 		}
		// 		cout<<endl;				
		// 	}
		// 	else
		// 		cout<<"has no conflict cluster"<<endl;

		// 	if(serialLocal != -1)
		// 	{			
		// 		cout<<"suvived conflict clusters:"<<endl;
		// 		int haveSuvivedConflict;
		// 		std::pair<int,int> toFindClusterPair;
		// 		for(int ITuncSET = 0; ITuncSET < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
		// 		{	
		// 			haveSuvivedConflict = find_ele.find(goodSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
		// 			if(haveSuvivedConflict != -1)
		// 			{
		// 				cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]<<" "<<endl;;

		// 				if(*cItSet > clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET])
		// 				{
		// 					pairToCheck.first  = *cItSet;
		// 					pairToCheck.second = clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET];
		// 				}

		// 				//it is a map that contains conflict clsuter pair and  its distance when add to cluter
		// 				toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET]);
		// 				if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
		// 					clusterizer.conflictClusterPairDistanceMap.end())
		// 				{
		// 					cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
		// 				}
		// 				else
		// 				{
		// 					cout<<"can not find conflict info in previous pair map"<<endl;
		// 					cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
		// 					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);

		// 					for(int ITuncSET_second = ITuncSET+1; ITuncSET_second < clusterizer.all_conflict_cluster[serialLocal].second.size(); ITuncSET_second++)
		// 					{
		// 						toFindClusterPair = std::pair<int , int > (*cItSet, clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]);
		// 						cout<<clusterizer.all_conflict_cluster[serialLocal].second[ITuncSET_second]<<" conflict cluster"<<endl;;

		// 						if(clusterizer.conflictClusterPairDistanceMap.find(toFindClusterPair) != 
		// 							clusterizer.conflictClusterPairDistanceMap.end())
		// 						{
		// 							cout<<"the distance when form this conflict is "<<clusterizer.conflictClusterPairDistanceMap[toFindClusterPair]<<endl;
		// 						}
		// 						else
		// 						{
		// 							cout<<"can not find conflict info in previous pair map"<<endl;
		// 							cout<<"the total number in the map is "<<clusterizer.conflictClusterPairDistanceMap.size()<<endl;
		// 							printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 						}	
		// 					}
		// 					exit(0);
		// 				}
		// 				cout<<"conflict exist after inter check"<<endl;
		// 				exit(0);
		// 			}
		// 		}
		// 		cout<<endl;				
		// 	}
		// 		//check if the uncertian cluster also exit in the survived cluster
		// }


		//intra and inter
		activeLoops.clear();
		gatherLinks(goodSet,activeLoops);
		const char *g2f_intra_inter=OptimizedResultIntraInter.data();
		gWrapper->saveOptimizedResult(g2f_intra_inter, activeLoops);
		int numofGoodS = activeLoops.size();


		// //intra
		// activeLoops.clear();
		// gatherLinks(suvivedTointerSet_backup,activeLoops);
		// const char *g2f_intra=OptimizedResultIntra.data();
		// gWrapper->saveOptimizedResult(g2f_intra, activeLoops);
		// int numoforiginal = activeLoops.size();
		// cout<<"loops in good set: "<<numofGoodS<<" and that in original: "<<numoforiginal<<endl;


		// //intra and conflict
		// activeLoops.clear();
		// gatherLinks(suvivedToInter,activeLoops);
		// const char *g2f_intra_conflict=OptimizedResultIntraConflict.data();
		// gWrapper->saveOptimizedResult(g2f_intra_conflict, activeLoops);
		// numofGoodS = activeLoops.size();


		// goodSet.clear();
		// interClusterConsistent_second( suvivedTointerSet_backup_sec, goodSet, rejectSet);

		// cout<<"cluster set for intra and conflict and inter: ";
		// 	cItSet	= goodSet.begin();
		// 	cEndSet = goodSet.end();
		// for( ; cItSet!=cEndSet; cItSet++)
		// {
		// 	cout<<*cItSet<<" ";
		// }
		// cout<<endl;

		// //intra and conflict and inter 
		// activeLoops.clear();
		// gatherLinks(goodSet,activeLoops);
		// const char *g2f_intra_conflict_inter=OptimizedResultIntraConflictInter.data();
		// gWrapper->saveOptimizedResult(g2f_intra_conflict_inter, activeLoops);
		// numoforiginal = activeLoops.size();
		// cout<<"loops in good set: "<<numofGoodS<<" and that in original: "<<numoforiginal<<endl;

		return 1;

		// std::set<int>  H_refine, H_refine_2;
		// multiClusterCheck(suvivedTointerSet, H_refine);//optimize and then delete all the clusters whose sum chi2 error is bigger the thereshold

		// activeLoops.clear();
		// gatherLinks(H_refine,activeLoops);
		// const char *g2f_=OptimizedResultWithRecheck.data();
		// gWrapper->saveOptimizedResult(g2f_, activeLoops);


		// H_refine.clear();
		// multiClusterCheck_onlyMinusOnePerTime(suvivedTointerSet, H_refine);//optimeze and delete only the biggest mean error one
		// activeLoops.clear();
		// gatherLinks(H_refine,activeLoops);
		// const char *g2f___=OptimizedResultWithReRecheck.data();
		// gWrapper->saveOptimizedResult(g2f___, activeLoops);	
		
		// cout<<"result after delete only the cluster who has the biggest chi2 error"<<endl;			
		// multiClusterCheck(H_refine, H_refine_2);

		// exit(0);

		//traversing the remined big clusters,check the uncertain set to see if its uncertain cluster also suvived, if suvived, combine the two cluster 
		//together and optimize, if both cluster meet the residual error limit, then they are consistent

		// std::vector<std::vector<std::vector<int>>> pyramid_group_cluster;
		// std::vector<int> lowLevel;
		// std::vector<std::vector<int> > middleLevel;

		// cIt 	= consistentClusters.begin(),
		// cEnd 	= consistentClusters.end();

		// for( ; cIt!=cEnd; cIt++)
		// {
		// 	lowLevel.clear();
		// 	lowLevel.push_back(*cIt);
		// 	middleLevel.push_back(lowLevel);
		// }
		// pyramid_group_cluster.push_back(middleLevel);

		// int conflictButPassTwoClusterCheck, TestClusterID, haveConflictInfo;

		// for(int level = 0; level < 1000; level++)
		// {
		// 	if(level == 0)//do two cluster check to determin the uncertain relationship between the suvived clusters 
		// 	{
		// 		middleLevel.clear();
		// 		for(int base = pyramid_group_cluster[0].size() - 1; base > 0; base--)
		// 		{
		// 			// int clusterIDtoFind = pyramid_group_cluster[level][base][0];
		// 			// int serialLocal = find_ele.find(clusterizer.all_uncertain_cluster, clusterIDtoFind);
		// 			for(int variable = base-1; variable >= 0; variable--)
		// 			{
		// 				// if(serialLocal != -1)
		// 				// {
		// 					// clusterIDtoFind = pyramid_group_cluster[level][variable][0];
		// 					// int Lssb = find_ele.find(clusterizer.all_uncertain_cluster[serialLocal].second, clusterIDtoFind);
		// 					// if(Lssb != -1)
		// 					{
		// 						hypotheses.clear();
		// 						hypotheses.insert(pyramid_group_cluster[0][base][0]);
		// 						hypotheses.insert(pyramid_group_cluster[0][variable][0]);
		// 						TestClusterID = pyramid_group_cluster[0][variable][0];

		// 						if(twoClusterCheck(hypotheses))
		// 						{
		// 							haveConflictInfo = find_ele.find(clusterizer.all_conflict_cluster,  pyramid_group_cluster[0][base][0]);
		// 							cout<<"pass find out if there is confict info to the base cluster"<<endl;
		// 							if(haveConflictInfo != -1)
		// 							{
		// 								conflictButPassTwoClusterCheck = find_ele.find(
		// 									clusterizer.all_conflict_cluster[haveConflictInfo].second, TestClusterID);
		// 								if(conflictButPassTwoClusterCheck != -1)
		// 								{
		// 									cout<<"two confilict cluster but all pass chi2 test in optimization process"<<endl;
		// 									printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 									exit(0);
		// 								}	
		// 							}

		// 							lowLevel.clear();
		// 							lowLevel.push_back(pyramid_group_cluster[level][base][0]);
		// 							lowLevel.push_back(pyramid_group_cluster[level][variable][0]);
		// 							// assert(lowLevel.size() == 2);
		// 							middleLevel.push_back(lowLevel);
		// 							cout<<"consistent clusters: "<<pyramid_group_cluster[level][base][0]<<" "<<
		// 								pyramid_group_cluster[level][variable][0]<<endl;
		// 						}
		// 						else
		// 						{
		// 							cout<<"confilict : "<<pyramid_group_cluster[level][base][0]<<" "<<
		// 								pyramid_group_cluster[level][variable][0]<<endl;
		// 						}

		// 					}
		// 				// }
		// 			}
		// 		}
				
		// 		pyramid_group_cluster.push_back(middleLevel);
		// 		// for(int ddd0 = 0; ddd0 < middleLevel.size(); ddd0++)
		// 		// {
		// 		// 	cout<<middleLevel[ddd0][0]<<" "<<middleLevel[ddd0][1]<<endl;
		// 		// 	if(middleLevel[ddd0].size() != 2)
		// 		// 		cout<<"elemetn not equal to 2"<<endl;
		// 		// 	// std::cout << "Press \'Return\' to end." << std::endl;
		// 		// 	// std::cin.get();
		// 		// }
		// 	}
		// 	else
		// 	{
		// 		middleLevel.clear();
		// 		for(int base = 0; base < pyramid_group_cluster[level].size(); base++) //group
		// 		{
		// 			// int serialLocal = find_ele.find(clusterizer.all_uncertain_cluster, pyramid_group_cluster[level][base]);


		// 			cout<<" "<<endl;
		// 			for(int variable = 0; variable < pyramid_group_cluster[0].size(); variable++)//cluster
		// 			{
		// 				cout<<" "<<endl;
		// 				if(find_ele.find(pyramid_group_cluster[level][base], pyramid_group_cluster[0][variable][0]) != -1)
		// 					continue;
		// 				bool add = 1;
		// 				lowLevel.clear();
		// 				lowLevel.push_back(pyramid_group_cluster[0][variable][0]);
		// 				cout<<"try to add: "<<pyramid_group_cluster[0][variable][0]<<endl;
		// 				cout<<"existing group: ";
		// 				for(int tower = 0; tower < pyramid_group_cluster[level][base].size(); tower++)
		// 				{
		// 					lowLevel.push_back(pyramid_group_cluster[level][base][tower]);
		// 					cout<<pyramid_group_cluster[level][base][tower]<<" ";
		// 				}
		// 				cout<<endl;
		// 				std::sort(lowLevel.begin(), lowLevel.end());
		// 				int ddlj = find_ele.find(middleLevel, lowLevel);
		// 				if(ddlj != -1)
		// 				{
		// 					// printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 					continue; 
		// 				}
		// 				// 

		// 				for(int tower = 0; tower < pyramid_group_cluster[level][base].size(); tower++)
		// 				{
		// 					std::pair<int, int> checkPair(pyramid_group_cluster[level][base][tower], pyramid_group_cluster[0][variable][0]);
		// 					int findCheckPair = find_ele.find(pyramid_group_cluster[1], checkPair);
		// 					if(findCheckPair == -1)
		// 					{
		// 						// cout<<"did not find cosistnet pair to support that it should be add to the group "<<endl;
		// 						cout<<"did not find consistent pair "<<checkPair.first<<" "<<
		// 						    checkPair.second<<endl;
		// 						add = 0;
		// 						break;
		// 					}
		// 					else
		// 					{
		// 						cout<<"find consistent pair "<<pyramid_group_cluster[1][findCheckPair][0]<<" "<<
		// 						    pyramid_group_cluster[1][findCheckPair][1]<<endl;
		// 					}
		// 				}
		// 				if(add == 1)
		// 				{
		// 					cout<<"consistent group: ";
		// 					for(int sd = 0; sd < lowLevel.size(); sd++)
		// 						cout<<lowLevel[sd]<<" ";
		// 					cout<<endl;
		// 					middleLevel.push_back(lowLevel);
		// 				}
		// 			}
		// 		}
		// 		if(middleLevel.size() == 0)
		// 		{
		// 			break;
		// 		}	
		// 		else	
		// 		{
		// 			pyramid_group_cluster.push_back(middleLevel);
		// 			for(int ddaa99 = 0; ddaa99 <  middleLevel.size(); ddaa99++)
		// 			{
		// 				for(int sd = 0; sd < middleLevel[ddaa99].size(); sd++)
		// 					cout<<middleLevel[ddaa99][sd]<<" ";
		// 				cout<<endl;	
		// 			}
		// 		}	
		// 	}

		// cout<<" "<<endl;
		// for(int dsflla093 = 0; dsflla093 < pyramid_group_cluster.size(); dsflla093++)
		// {
		// 	cout<<"level["<< dsflla093<<"] has "<<pyramid_group_cluster[dsflla093].size()<<" elements"<<endl;
		// }

		// std::cout << "Press \'Return\' to end." << std::endl;
		// std::cin.get();
		// }
		// cout<<" "<<endl;
		// for(int dsflla093 = 0; dsflla093 < pyramid_group_cluster.size(); dsflla093++)
		// {
		// 	cout<<"level["<< dsflla093<<"] has "<<pyramid_group_cluster[dsflla093].size()<<" elements"<<endl;
		// }

		// cout<<" "<<endl;
		// cout<<"pyramid_group_cluster.back().size(): "<<pyramid_group_cluster.back().size()<<endl;
		// for(int lkhh3 =0; lkhh3 < pyramid_group_cluster.back().size(); lkhh3++)
		// {
		// 	int mutiOptimize = multiClusterCheck(pyramid_group_cluster.back()[lkhh3]);
		// }
		
		// exit(0);


		cItSet 	= hypotheses_largeCluster.begin();
		cEndSet 	= hypotheses_largeCluster.end();
		std::vector<std::set<int> > all_regain_info;
		all_regain_info.clear();
		cout<<" "<<endl;
		cout<<"************************************"<<endl;
		cout<<"start to merge uncertian clusters: "<<endl;
		for( ; cItSet!=cEndSet; cItSet++)//traverse all survived cluster
		{
			cout<<"*********************************************8 "<<endl;
			cout<<"current base cluster is"<<*cItSet<<endl;
			//make sure the suvived clsuter has uncertain cluster, *cItSet is the serial number of the test cluster
			int serialLocal = find_ele.find(clusterizer.all_uncertain_cluster, *cItSet);
			if(serialLocal != -1)//if the current cluster has some uncertain clusters
			{
				cout<<"this base cluster has uncertain clusters"<<endl;
				std::set<int> mergeSEG;
				mergeSEG.clear();
				for(int ITuncSET = 0; ITuncSET < clusterizer.all_uncertain_cluster[serialLocal].second.size(); ITuncSET++)//tranverse each element in the uncertain set of the current cluster
				{	
					// if(consistentClusters.find(clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]) != consistentClusters.end());//check if *ITuncSET has suvived
					int NoRepeatSerialLocal = find_ele.find(consistentClusters, clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]);
					if(NoRepeatSerialLocal != -1)
					{
						cout<<"uncertain cluster that suvived of the base cluster: "<<clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]<<endl;
						cout<<"the founded cluster is: "<<consistentClusters[NoRepeatSerialLocal]<<endl;
						if(consistentClusters[NoRepeatSerialLocal] != 
							clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET])
						{
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
						}
						cout<<" "<<endl;
						// if(*(consistentClusters.find(*ITuncSET)) != *ITuncSET)
						// 	continue;
						hypotheses.clear();
						hypotheses.insert(clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]);
						hypotheses.insert(*cItSet);
						if(*cItSet == clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET])
							continue;
						cout<<"two cluster check involve clusters: "<<clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]<<" and "<<*cItSet<<endl;
						if(twoClusterCheck(hypotheses))
						{
							mergeSEG.insert(clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET]);
							mergeSEG.insert(*cItSet);
							cout<<"cluster "<< (clusterizer.all_uncertain_cluster[serialLocal].second[ITuncSET])<<"get a good uncertain cluster: "<<*cItSet<<endl;
						}
						else
						{
							cout<<"its a bad cluster pair"<<endl;
						}
					}
				}

				int its;
				std::vector<std::pair<int, int>> conflictClustersGainedFromUncertain, consistentClustersGainedFromUncertain;

				//this merge process should consider the 
				if(mergeSEG.size()!= 0)//
				{
					bool bit =0;
					std::set<int> & foradd = mergeSEG;
					std::set<int>::iterator it = mergeSEG.begin(), iend = mergeSEG.end();
					for(; it != iend; it++)
					{

						// std::set<std::set<int> >::iterator its = all_regain_info.begin(), iends = all_regain_info.end();
						for(its = 0; its < all_regain_info.size(); its++)
						{
							// if((*its).find(*it) != (*its).end())
							if(all_regain_info[its].find(*it) != all_regain_info[its].end()) 
							{
								bit    = 1;
								// foradd = *its;
								// cout<<"type of *its: "<<typeof(*its)<<endl;
								all_regain_info[its].insert(mergeSEG.begin(), mergeSEG.end());
								// all_regain_info[*its.first].insert(mergeSEG.begin(), mergeSEG.end());
								break;
							}
						}
						if(bit == 1)
							break;
					}

					it = mergeSEG.begin(), iend = mergeSEG.end();
					cout<<"mergeSEG: ";
					for(; it != iend; it++)
					{
						cout<<*(it)<<" ";
					}
					cout<<endl;
					if(bit == 1)
					{
						cout<<"current group:";
						std::set<int >::iterator iteratorPrint = all_regain_info[its].begin(), itpringend = all_regain_info[its].end();
						for(; iteratorPrint != itpringend; iteratorPrint++)
						{
							cout<<*(iteratorPrint)<<" "<<endl;
						}
						// foradd.insert(mergeSEG.begin(), mergeSEG.end());
						cout<<endl;
					}
					else
					{
						all_regain_info.push_back(mergeSEG);
					}
					cout<<"cluster "<< *cItSet<<" has consistent cluster after uncertain check"<<endl;
					// all_regain_info.insert(mergeSEG);
				}
			}
			else
				cout<<"cluster "<< *cItSet<<" has no uncertain cluster"<<endl;
		}

		// exit(0);
		std::set<int> totalMerge;
		// cout<<"all_regain_info size: "<<all_regain_info.size()<<endl;
		for(int mit =0; mit < all_regain_info.size(); mit++)
		{
			// cout<<" regain "<<*mit<<endl;
			totalMerge.insert(all_regain_info[mit].begin(), all_regain_info[mit].end());
		}
		// cout<<"totalMerge size: "<<totalMerge.size()<<endl;

		std::vector<std::vector<int> > finalVertor;
		std::set<int> proposed_mergeCluseter;
		// std::set<int>::iterator heit = consistentClusters.begin(), heenhd = consistentClusters.end();
		std::vector<int>::iterator heit = consistentClusters.begin(), heenhd = consistentClusters.end();
		for(; heit != heenhd; heit++)//for all suvived clusters
		{

			if(proposed_mergeCluseter.find(*heit) != proposed_mergeCluseter.end())// to be improved
				continue;
			std::vector<int> single_final_cluster;
			if(totalMerge.find(*heit) != totalMerge.end())// if this cluster has to be merged cluster
			{
				// cout<<" if"<<endl;
				
				bool bit =0;
				for(int mit =0; mit < all_regain_info.size(); mit++)//
				{
					if(all_regain_info[mit].find(*heit) != all_regain_info[mit].end())//find the set of the merge cluster
					{
						std::set<int>::iterator innit = all_regain_info[mit].begin(), innend =all_regain_info[mit].end();
						for(; innit != innend; innit++)
						{
							single_final_cluster.push_back(*innit);
							proposed_mergeCluseter.insert(*innit);
						}
						bit = 1;
						// cout<<"if"<<endl;
						finalVertor.push_back(single_final_cluster);
						break;
					}
				}
				if(bit == 0)
				{
					cout<<"should find merge cluster, but not"<<endl;
					exit(0);
				}
			}
			else
			{
				// cout<<"else"<<endl;
				single_final_cluster.push_back(*heit);
				finalVertor.push_back(single_final_cluster);
			}
		}
		// cout<<" finalVertor size: "<<finalVertor.size()<<endl;
	
		std::vector<std::pair<int, std::vector<int> > > gropuANDloopNum;
		std::pair<int, std::vector<int> > fiele;
		bool bit4nonsingleElementCluster = 0;
		for(int i = 0; i < finalVertor.size(); i++)
		{
			int totn = 0;//the total number of the group
			for(int j = 0; j < finalVertor[i].size(); j++)
			{
				totn = totn+clusterizer.getClusterByID(finalVertor[i][j]).size();
			}
			fiele.first  = totn;
			if (totn > 1)
				bit4nonsingleElementCluster = 1;
			fiele.second = finalVertor[i];
			gropuANDloopNum.push_back(fiele);
		}

		// if(finalVertor.size() == 0)
		// {

		// }
		
		cout<<" gropuANDloopNum size: "<<gropuANDloopNum.size()<<endl;
		cout<<" finalVertor size: "<<finalVertor.size()<<endl;

		//if no loop suvived, save odometry info only and return 0
		if(finalVertor.size() == 0)
		{
			activeLoops.clear();
			// gWrapper->optimize(activeLoops,nIterations, odoEdgeRelateLC_Error, odoEdgeError);
			//save the result file before doube check 
			// gWrapper->optimizer->save(OptimizedResultWithoutRecheck);
			const char *g2f=OptimizedResultWithoutRecheck.data();
			gWrapper->saveOptimizedResult(g2f, activeLoops);
			cout<<"no loops suvived"<<endl;
			if(clusterizer.set4GTL.size()==0)
			{
				cout<<"rignt, there is indeed not ground truth loop closure"<<endl;
			}
			else
				cout<<"wrong , as it failed to pick out the true LC"<<endl;
			sleep(3);
			return 0;
		}
		// sleep(3);
		// cout<<"debug 3"<<endl;
		if(bit4nonsingleElementCluster == 1)
			std::sort(gropuANDloopNum.begin(), gropuANDloopNum.end(), cmp_vertor_int_pair);
		else
		{
			for(int i = 0; i < gropuANDloopNum.size();  i++)
			{
				cout<<"gropuANDloopNum[i].second.size: "<<gropuANDloopNum[i].second.size()<<endl;
				if(gropuANDloopNum[i].second.size() != 1)
				{
					cout<<"the size should be 1"<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				gropuANDloopNum[i].first = floor((overallChi2SET4singleElementCluster[gropuANDloopNum[i].second[0]])*(-10));
			}
			std::sort(gropuANDloopNum.begin(), gropuANDloopNum.end(), cmp_vertor_int_pair);
		}
		// activeLoops.clear();
		// gatherLinks(hypotheses_largeCluster,activeLoops);
		// linkErrors =gWrapper->optimize(activeLoops,nIterations);
		// cout<<" "<<endl;
		// cout<<"re optimise on clusters who has menber more than one: "<<endl;
		// std::set<int>::iterator bigCluster = hypotheses_largeCluster.begin(), bCend = hypotheses_largeCluster.end();
		// for(;bigCluster !=  bCend; bigCluster++)
		// 	cout<<*bigCluster<<endl;
		// for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
		// 		it!=end;
		// 		it++)
		// {
		// 	if(it->first.first < 0) continue;
		// 	// if(clusterizer.getClusterID(it->first) == 181)
		// 		cout<<"loop "<<it->first.first << " "<< it->first.second <<" in " <<"cluster "<<clusterizer.getClusterID(it->first)<<
		// 			"has residual error "<<it->second<<endl;
		// 	if( it->second < utils::chi2(edgeDimension))
		// 	{
		// 		hypotheses.insert(clusterizer.getClusterID(it->first));
		// 		cout<<clusterizer.getClusterID(it->first)<<endl;
		// 		done = false;
		// 	}
		// }

		// exit(0);

		hypotheses.clear();		
		std::cout<<" GoodSet :";
		// for( IntSet::iterator it= goodSet.begin(), end = goodSet.end(); it!=end; it++)
		cout<<"ele number in the firs group :"<<gropuANDloopNum[0].second.size()<<endl;

		cout<<" "<<endl;
	 cout<<"the biggest group menber cluster: ";
		for( std::vector<int>::iterator it= gropuANDloopNum[0].second.begin(), end = gropuANDloopNum[0].second.end(); it!=end; it++)
		{
	
			std::cout<<(*it)<<" ";
			hypotheses.insert(*it);
		}
		std::cout<<std::endl;
		//print the [0]
		if(gropuANDloopNum.size() > 1)
		{
			cout<<" "<<endl;
			 cout<<"the second biggest group menber cluster: ";
			for( std::vector<int>::iterator it= gropuANDloopNum[1].second.begin(), end = gropuANDloopNum[1].second.end(); it!=end; it++)
			{
		
				std::cout<<(*it)<<" ";
				// hypotheses.insert(*it);
			}
			std::cout<<std::endl;
		}

		//ACCURACY : check if all suvived loops are good one
		activeLoops.clear();
		gatherLinks(hypotheses,activeLoops);

		cout<<"activeLoops num which should equal to that in the suvived: "<<activeLoops.size()<<endl;

		//final optimize on the all loops
		IntPairDoubleMap linkErrors;
		cout<<" "<<endl;
		cout<<"cluster ID of qualified loops : "<<endl;
		std::vector<std::pair<std::pair<int,int>, int> >final_bad_loop;
		std::pair<std::pair<int,int>, int> ele_final_bad_loop;
		int biggestERRIRclusterSerial = -1;

		if(gropuANDloopNum[0].second.size() > 1)
		{
			linkErrors =gWrapper->optimize(activeLoops,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			int error =0;
			for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
					it!=end;
					it++)
			{
				if(it->first.first < 0) continue;
		
				cout<<"loop "<<it->first.first<<" "<<it->first.second<<" in cluster "<< clusterizer.getClusterID(it->first)<<" has error "<<it->second<<endl;
				// if( it->second > utils::chi2(edgeDimension))		
				if( it->second > twoNineBelief)
				{
					// cout<<clusterizer.getClusterID(it->first)<<endl;
					ele_final_bad_loop.first  = it->first;
					ele_final_bad_loop.second = clusterizer.getClusterID(it->first);
					final_bad_loop.push_back(ele_final_bad_loop);
					cout<<" "<<endl;
					if (error < it->second)
					{
						error = it->second;
						biggestERRIRclusterSerial = ele_final_bad_loop.second;
					}
					// done = false;
				}
			}
		}

		//if there are some  bad loop from a cluster whose memner loop number is only one, delete it and reoptimize
		bool proposed = 0;
		for(int i = 0; i < final_bad_loop.size(); i++)
		{
			if(clusterizer.getClusterByID(final_bad_loop[i].second).size() == 1)
			{
				if(final_bad_loop[i].second == biggestERRIRclusterSerial)
				{
					hypotheses.erase(final_bad_loop[i].second);	

					proposed = 1;
				}

			}
			// else
			// {
			// 	doubtLoopVector.push_back(final_bad_loop[i].first);
			// 	clusterIDofDoubtLoops.push_back(final_bad_loop[i].second);
			// }
		}

		// if(proposed == 0 and final_bad_loop.size() != 0)
		// {
		// 	cout<<"the biggest error loops is not in a single element cluster"<<endl;
		// 	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// 	exit(0);
		// }
		//if 
		// final_bad_loop.clear();
		if(hypotheses.size() > 1 and final_bad_loop.size() != 0)
		{
			final_bad_loop.clear();
	 		activeLoops.clear();
			gatherLinks(hypotheses,activeLoops);
			linkErrors =gWrapper->optimize(activeLoops,small_numIterations, odoEdgeRelateLC_Error, odoEdgeError);
			for ( IntPairDoubleMap::iterator it = linkErrors.begin(), end = linkErrors.end();
					it!=end;
					it++)
			{
				if(it->first.first < 0) continue;
		
				cout<<"loop "<<it->first.first<<" "<<it->first.second<<" in cluster "<< clusterizer.getClusterID(it->first)<<" has error "<<it->second<<endl;
				// if( it->second > 50)
				// {
				// 	// cout<<clusterizer.getClusterID(it->first)<<endl;
				// 	ele_final_bad_loop.first  = it->first;
				// 	ele_final_bad_loop.second = clusterizer.getClusterID(it->first);
				// 	final_bad_loop.push_back(ele_final_bad_loop);
				// 	cout<<" "<<endl;
				// 	// done = false;
				// }
				// else 
				if( it->second > utils::chi2(edgeDimension))
				{
					// cout<<clusterizer.getClusterID(it->first)<<endl;
					ele_final_bad_loop.first  = it->first;
					ele_final_bad_loop.second = clusterizer.getClusterID(it->first);
					final_bad_loop.push_back(ele_final_bad_loop);
					cout<<" "<<endl;
					// done = false;
				}
			}
		}


		//if the bad loop detected from the last step comes from a cluster whose memner loop number is only one,
		//then delete the cluster

		for(int i = 0; i < final_bad_loop.size(); i++)
		{
			if(clusterizer.getClusterByID(final_bad_loop[i].second).size() == 1)
			{
				hypotheses.erase(final_bad_loop[i].second);
			}
			else
			{
				cout<<"the final bad loop "<<final_bad_loop[i].first.first<<" "<<final_bad_loop[i].first.second
					<<" in cluster "<< clusterizer.getClusterID(final_bad_loop[i].first)<<" comes form a cluster whose menber loops are more than one"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);

				doubtLoopVector.push_back(final_bad_loop[i].first);
				clusterIDofDoubtLoops.push_back(final_bad_loop[i].second);

				clusterizer.setClusterID(final_bad_loop[i].first,ID_IGNORE);
				// exit(0);
			}
		}

		//after delete the bad cluster, re get all suvived loops
		activeLoops.clear();
		gatherLinks(hypotheses,activeLoops);
		// gWrapper->optimize(activeLoops,nIterations, odoEdgeRelateLC_Error, odoEdgeError);
		//save the result file before doube check 
		// gWrapper->optimizer->save(OptimizedResultWithoutRecheck);
		const char *g2f=OptimizedResultWithoutRecheck.data();
		gWrapper->saveOptimizedResult(g2f, activeLoops);

		// //save the good loops into a g2o file
		// ofstream fileStreamr; 
		// fileStreamr.open(final_suvived_loops_fileName,ios::trunc);
		// std::pair<g2o::SE2, Matrix3d> tSave;
		// std::pair<int,int> ty;
		// int  originalLoopNum;
		// std::vector<std::pair<double, std::pair<int, int> > >  chiStstis4eachLoop;
		// std::vector<std::pair<int, int> >  doubtLoopVector;
		// std::vector<int> clusterIDofDoubtLoops;

		// for(size_t i=0; i< clusterizer.clusterCount(); i++)
		// {
		// 	std::cout<<i<<" "; std::cout.flush();
		// 	std::vector<double>  chiStatisVector;
		// 	std::vector<std::pair<int, int> >  doubtLoopNumberVector;
		// 	if(intraClusterConsistent(i, chiStatisVector, doubtLoopNumberVector, originalLoopNum, chiStstis4eachLoop))
		// 	{
		// 		for(int doubt = 0; doubt < doubtLoopNumberVector.size(); doubt++)
		// 		{

		// 			doubtLoopVector.push_back(doubtLoopNumberVector[doubt]);
		// 			clusterIDofDoubtLoops.push_back(i);
		// 		}
		// 		//save loops information in i to a txt file
		// 		for(int cbb = 0; cbb < clusterizer._clustersFound[i].positionserial.size(); cbb++)
		// 		{
		// 			ty.first = clusterizer._clustersFound[i].positionserial[cbb][0];
		// 			ty.second = clusterizer._clustersFound[i].positionserial[cbb][3];

		// 			tSave = clusterizer.LP_Trans_Covar_Map[ty];
		// 			Matrix3d fsave;
		// 			fsave = tSave.second;

		// 			fileStreamr<<"EDGE_SE2 "<<ty.first<<" "<<ty.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
		// 				<<" "<<fsave.inverse()(0,0)<<" "<<fsave.inverse()(0,1)<<" "<<fsave.inverse()(0,2)<<" "<<
		// 				fsave.inverse()(1,1)<<" "<<fsave.inverse()(1,2)<<" "<<fsave.inverse()(2,2)<<" "<<"\n";
		// 		}
		// 		consistentClusters.insert(i);
		// 		loopNUmbersInConsistentClusters.push_back(std::pair<int, int> (i, originalLoopNum));
		// 	}
		// }
		// 		cout<<"consistentClusters size: "<<consistentClusters.size()<<endl;
		// 		sleep(2);
		// fileStreamr.close();


		int  forRecall=0;// number of the true loops that has been successfully picked out
		int goodLN = 0;  // number of good loops in all the picked out loops
		IntPairSet::iterator iit = activeLoops.begin(), iiend = activeLoops.end();
		int allT = 0;

		//accuracy
		if(clusterizer.NumOfRealLoops != -1)
		{
			for(; iit != iiend; iit++)
			{
				if(clusterizer.set4GTL.find(*iit) == clusterizer.set4GTL.end())
				{
					allT = allT+1;
					cout<<"suvived loop "<<(*iit).first <<" "<<(*iit).second<< " is a flase one"<<endl;
				}
				else
					goodLN = goodLN+1;
			}
			if(allT != 0)
				cout<<" there are "<< allT<<" bad loops"<<endl;
			else
				cout<<"all loops suvived are true "<<endl;
		}

		//recall: find out the true loop that was abandened
		if(clusterizer.NumOfRealLoops != -1)
		{		
			// if(clusterizer.set4GTL.size() > activeLoops.size())
			// {
				iit = clusterizer.set4GTL.begin(), iiend = clusterizer.set4GTL.end();
			
				for(; iit != iiend; iit++)
				{
					if(activeLoops.find(*iit) == activeLoops.end())
					{
						cout<<"true loop "<<(*iit).first <<" "<<(*iit).second<< " is abandened"<<endl;
					}
					else
						forRecall = forRecall+1;

				}
			// }
		}

		//if there exists an ground truth date, then print out the relevent information
		cout<<" "<<endl;
		if(clusterizer.NumOfRealLoops != -1)
		{
			cout<<" ******** information without recheck doubt loops ******** "<<endl;
			cout<<"the number of     suvived good loops    is: "<<activeLoops.size()<<std::endl;
			cout<<" "<<endl;
			cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the                  accuracy rate      is: "<<
				(goodLN)*1.0/(activeLoops.size())<<" "<<endl;
			cout<<"goodLN "<< goodLN <<endl;
			cout<<" "<<endl;
			cout<<"the                  recall  rate       is: "<< 
				double(forRecall*1.0)/(clusterizer.set4GTL.size())<<" "<<endl;
				cout<<"forRecall: "<<forRecall<<endl;
			cout<<" ********************************************************* "<<endl;
			cout<<" "<<endl;
		}

		//recheck the doubt loops
		std::vector<int> trueBitOfReacceptDoubtLoop;
		if (doubtLoopVector.size() > 0)
		{
			if(clusterIDofDoubtLoops.size() != doubtLoopVector.size())
			{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}
			
			for(int doubt = 0; doubt < doubtLoopVector.size(); doubt++)
			{
				cout<<"doubt chenck "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second;
				std::vector< std::pair<double, std::pair<std::pair<int, int>, int> > > & middCheck = 
					clusterizer._clustersFound[clusterIDofDoubtLoops[doubt]].chiStatisWhenAdd2Cluster;;
				// middCheck = clusterizer._clustersFound[clusterIDofDoubtLoops[doubt]].chiStatisWhenAdd2Cluster;
				bool find = 0;
				int j=0;
				for(j =0; j< middCheck.size(); j++)
				{
					if(middCheck[j].second.first == doubtLoopVector[doubt])
					{
						find = 1;
						cout<<" has transdis "<<middCheck[j].first<<endl;
						if(middCheck[j].first < chi2_continuous(edgeDimension, 0.6))
						{
							reAccept.push_back(doubtLoopVector[doubt]);
							//reput the loop into the cluster
							clusterizer.setClusterID(doubtLoopVector[doubt],clusterIDofDoubtLoops[doubt]);
							//check the reaccept are true or not

							if(clusterizer.set4GTL.find(doubtLoopVector[doubt]) == clusterizer.set4GTL.end())
								trueBitOfReacceptDoubtLoop.push_back(0);
							else
								trueBitOfReacceptDoubtLoop.push_back(1);

						}
						else
						{
							cout<<" this loop "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second<<" has been deleted"<<endl;
							clusterizer.setClusterID(doubtLoopVector[doubt],ID_IGNORE);
						}
						break;
					}
				}
				if(find == 0)
				{
					cout<<"the loop in doubt should be find out in the previous cluster, however, falied!"<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}

				cout<<"doubt loop "<<doubt<<" is "<<doubtLoopVector[doubt].first<<" "<<doubtLoopVector[doubt].second<<
					" ,  the chi statis is "<<middCheck[j].first<< " when it is added to the cluster "<<
					clusterIDofDoubtLoops[doubt]<<endl;
			}	
		}
		else
			cout<<"no doube loop need to be checked"<<endl;

		//print out the information of reaccept loops,
		int addTrue = 0, addFalse = 0;
		for(int reacce = 0; reacce < reAccept.size(); reacce++)
		{
			cout<<"reaccept doubt loop "<<reAccept[reacce].first<<" "<<reAccept[reacce].second;

			if(trueBitOfReacceptDoubtLoop[reacce])
			{
				addTrue =addTrue +1;
				cout<<" is true"<<endl;
			}
			else
			{
				addFalse =addFalse +1;
				cout<<" is false"<<endl;
			}
		}


		//if there exists an ground truth date, then print out the relevent information
		cout<<" "<<endl;
		if((clusterizer.NumOfRealLoops != -1) and (reAccept.size()>0))
		{
			cout<<" ***************  after recheck doubt loops *************** "<<endl;
			cout<<"the number of     suvived good loops    is: "<<activeLoops.size()+addTrue+addFalse<<std::endl;
			cout<<" "<<endl;
			cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the                  accuracy rate      is: "<<
				(goodLN+addTrue)*1.0/(activeLoops.size()+addTrue+addFalse)<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the                  recall  rate       is: "<< 
				(forRecall+addTrue)*1.0/(clusterizer.set4GTL.size())<<" "<<endl;
			cout<<" ********************************************************* "<<endl;
			cout<<" "<<endl;
		}
		else
			cout<<"no change happen after recheck the doubt loops "<<endl;

		//re get the loops in activeLoops after boubt check
		activeLoops.clear();
		gatherLinks(hypotheses,activeLoops);
				//save the result file before doube check 
		// gWrapper->optimize(activeLoops,nIterations, odoEdgeRelateLC_Error, odoEdgeError);

		if((reAccept.size()>0))
		{
			const char *g2fre=OptimizedResultWithRecheck.data();
			gWrapper->saveOptimizedResult(g2fre, activeLoops);
		}
		//save the cluster info after doubt check
		// ofstream fileStreamr; 
		fileStreamr.open("loops_info_of_all_suvived_clusters_after_intracheck_after_doubtcheck.txt",ios::trunc);

		cItSet 	= hypotheses.begin(),
		cEndSet 	= hypotheses.end();
		cout<<"remained clusters in hypotheses: "<<endl;
		for(; cItSet != cEndSet; cItSet++)
			cout<<*cItSet<<" ";
		cout<<endl;
		cout<<" "<<endl;
		cout<<"clusters in consistentClusters"<<endl;
		cIt 	= consistentClusters.begin(),
		cEnd 	= consistentClusters.end();
		for(; cIt != cEnd; cIt++)
		{
			std::cout<<*cIt<<" "; std::cout.flush();
			std::vector<double>  chiStatisVector;
			std::vector<std::pair<int, int> >  doubtLoopNumberVector;

				fileStreamr<<*cIt<<"\n";
				//save loops information in i to a txt file
				IntPairSet::iterator it = clusterizer.getClusterByID(*cIt).begin(), itend = clusterizer.getClusterByID(*cIt).end();
				for(; it != itend; it++)
				{
					int first = (*it).first, toNode = (*it).second;
					int j;
					for(j = 0; j<clusterizer._clustersFound[*cIt].positionserial.size();j++)
					{
						if( first == clusterizer._clustersFound[*cIt].positionserial[j][0] and
							toNode == clusterizer._clustersFound[*cIt].positionserial[j][3])
						{
							ty.first = clusterizer._clustersFound[*cIt].positionserial[j][0];
							ty.second = clusterizer._clustersFound[*cIt].positionserial[j][3];

							tSave = clusterizer.LP_Trans_Covar_Map[ty];
							Matrix3d fsave;
							fsave = tSave.second;

							fileStreamr<<"EDGE_SE2 "<<ty.first<<" "<<ty.second<<" "<<tSave.first[0]<<" "<<tSave.first[1]<<" "<<tSave.first[2]
								<<" "<<fsave.inverse()(0,0)<<" "<<fsave.inverse()(0,1)<<" "<<fsave.inverse()(0,2)<<" "<<
								fsave.inverse()(1,1)<<" "<<fsave.inverse()(1,2)<<" "<<fsave.inverse()(2,2)<<" "<<"\n";	
							break;
						}
						if (j == clusterizer._clustersFound[*cIt].positionserial.size()-1)
						{
							cout<<"can not find the postion info from the cluster for the suvived loop "<<first<<" "<<toNode <<endl;
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
						}
	
					}
				}				
			
		}

		fileStreamr.close();
		cout<<" "<<endl;
		if((clusterizer.NumOfRealLoops != -1) and (reAccept.size()>0))
		{

			cout<<" ***************  after recheck doubt loops, to see if the reaccept loops have been reput into cluster *************** "<<endl;
			cout<<"the number of     suvived good loops    is: "<<activeLoops.size()<<std::endl;
			cout<<" "<<endl;
			cout<<"the number of         true  loops       is: "<< clusterizer.NumOfRealLoops<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the number of      loops in testfile    is: "<< clusterizer.LP_Trans_Covar_Map.size()<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the                  accuracy rate      is: "<<
				(goodLN+addTrue)*1.0/(activeLoops.size())<<" "<<endl;
			cout<<" "<<endl;
			cout<<"the                  recall  rate       is: "<< 
				(forRecall+addTrue)*1.0/(clusterizer.set4GTL.size())<<" "<<endl;
			cout<<" ********************************************************* "<<endl;
			cout<<" "<<endl;
		}
		else
			cout<<"no change happen after recheck the doubt loops "<<endl;
		_goodSet.insert(gropuANDloopNum[0].second.begin(), gropuANDloopNum[0].second.end());
		// sleep(3);
		return true;
	}

	bool write(const char* filename)
	{
		IntPairSet correctLinks ;
		gatherLinks(_goodSet, correctLinks);
		gWrapper->write(filename,correctLinks);

		return true;
	}

	bool write_resolved_result(const char* filename, const char* resultFileName)
	{
		IntPairSet correctLinks ;
		gWrapper->saveOptimizedResult(resultFileName);
		gatherLinks(_goodSet, correctLinks);
		gWrapper->write(filename,correctLinks);
		return true;
	}

	bool removeIncorrectLoops()
	{
		size_t count = clusterizer.clusterCount();

		IntSet rejected;

		for(size_t i = 0 ; i<count ; i++)
		{
			if(_goodSet.find(i)==_goodSet.end())
			{
				rejected.insert(i);
			}
		}
		rejected.insert(ID_IGNORE);

		IntPairSet rejectedLinks;

		gatherLinks(rejected,rejectedLinks);
		gWrapper->removeEdges(rejectedLinks);

		// Clear up the clusters from Clusterizer

		for( IntSet::const_iterator it = rejected.begin(), end = rejected.end(); it!=end; it++)
		{
			clusterizer.deleteCluster(*it);
		}

		return true;
	}

};




#endif /* RRR_HPP_ */
