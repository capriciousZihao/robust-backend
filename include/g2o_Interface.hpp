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


#ifndef BACKEND_G2O_HPP_
#define BACKEND_G2O_HPP_

#include "backEndInterface.hpp"

#include "g2o/types/slam2d/vertex_se2.h"

#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/factory.h"
#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
#include <Eigen/Geometry>
#include <numeric>
#include <stdio.h>

#include<iostream>

using namespace std;

typedef g2o::BlockSolver< g2o::BlockSolverTraits<-1, -1> >  		SlamBlockSolver;
typedef g2o::LinearSolverCSparse<SlamBlockSolver::PoseMatrixType> 	SlamLinearSolver;

template<typename VertexType, typename EdgeType>
class G2O_Interface : public BaseBackend
{
	typedef
			std::map< IntPair, g2o::HyperGraph::Edge* >
			IntPairToEdgePtrMap;

	typedef
			std::vector< g2o::HyperGraph::Edge* >
			EdgePtrVector;

	typedef g2o::SparseOptimizer OptimizerType;

	EdgePtrVector 					odomteryEdges;
	IntPairToEdgePtrMap				loopclosureEdges;

	
	g2o::OptimizableGraph::EdgeSet 	activeEdges, odoEdgeSegment;
	g2o::OptimizableGraph::VertexSet 	activeVertexes;
	// Transform<double, 3, 1>         try_to_get_vertex_position;
	Eigen::Matrix<double, 3, 3>     M33d;


	bool initialized;

public:
	g2o::SparseOptimizer*			optimizer;
	int odoSize;
	G2O_Interface()
	{
		optimizer = NULL;
		initialized = false;
	}

	bool setOptimizer(void* opt)
	{
		if(optimizer==NULL)
		{
			optimizer = (OptimizerType*)opt;
			initialized = true;
			store();
		}
		else
		{
			std::cerr<<"Already existing optimizer?"<<std::endl;
			return false;
		}
		return true;
	}


	bool getLoopClosures(IntPairSet& loops)
	{
		if(optimizer == NULL)
		{
			std::cerr<<"Please read in a g2o file or pass the pointer to an existing optimizer before calling getLoopClosures()"<<std::endl;
			return false;
		}


		g2o::OptimizableGraph::EdgeSet::iterator
		eIt = optimizer->edges().begin(),
		eEnd = optimizer->edges().end();

		for( ; eIt!=eEnd ; eIt++)
		{
			int e1 = (*eIt)->vertices()[0]->id();
			int e2 = (*eIt)->vertices()[1]->id();
			if(std::abs(e1-e2) > 1)
			{
				// if(e1 == 2480)
				// {
				// 	cout<<"find the desired node!"<<endl;
				// 	exit(0);
				// }
				loops.insert(IntPair(e1,e2));
				loopclosureEdges[IntPair(e1,e2)] = *eIt;

			}
			else
			{
				odomteryEdges.push_back(*eIt);
			}
		}
		odoSize = odomteryEdges.size();
		// cout<<"loop size: "<<loops.size()<<endl;
		// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
		// exit(0);

		for(int i=0; i<odomteryEdges.size(); i++)
		{
			if(odomteryEdges[i]->vertex(0)->id() != i or odomteryEdges[i]->vertex(1)->id() != (i+1))
			{
				cout<<"check if it is odometry edge fail"<<endl;
				cout<<"i is "<<i<<" but first vertex id is "<< 
				 	odomteryEdges[i]->vertex(0)->id()<<" and second vertex id is "<<odomteryEdges[i]->vertex(1)->id()<<endl;
				exit(0);
			}
		}

		return true;
	}

	virtual ~G2O_Interface(){}

	int vertexCount() { return optimizer->vertices().size(); };
	int edgeCount() { return optimizer->edges().size(); };
	


	bool store() // store the current graph state
	{
		g2o::OptimizableGraph::VertexIDMap::iterator
			vIt = optimizer->vertices().begin(),
			vEnd = optimizer->vertices().end();

		for(; vIt!=vEnd ; vIt++ )	static_cast< VertexType* >(vIt->second)->push();

		return true;
	}

	bool restore() // restore from backup the previous estimate of the graph
	{
		g2o::OptimizableGraph::VertexIDMap::iterator
			vIt = optimizer->vertices().begin(),
			vEnd = optimizer->vertices().end();

		for(; vIt!=vEnd ; vIt++ )	static_cast< VertexType* >(vIt->second)->pop();

		// HACK : store again!
		store();
		return true;

	}

	bool read(
			const char* filename
			)
	{
		if(optimizer != NULL)
		{
			std::cerr<<"An allocated optimizer already exists "<<std::endl;
			return false;
		}
		optimizer = new g2o::SparseOptimizer;
		//SlamLinearSolver* linearSolver = new SlamLinearSolver();
		std::unique_ptr<SlamLinearSolver> linearSolver (new SlamLinearSolver());
		linearSolver->setBlockOrdering(false);
		//SlamBlockSolver* blockSolver = new SlamBlockSolver(linearSolver);
		std::unique_ptr<SlamBlockSolver> blockSolver ( new SlamBlockSolver (std::move(linearSolver)));
		//g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(blockSolver);
		g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(std::move(blockSolver));
		optimizer->setAlgorithm(solverGauss);


		if(!optimizer->load(filename)){
			std::cerr<<"Can't find file to read : "<<filename<<std::endl;;
			return false;
		}
		store();
		initialized = true;
		return true;
	}

	// Optimize on a given set of loop closure links
	// @return: The error of each link + overall error of the active part of the graph as the last element of the vector

	// IntPairDoubleMap optimize(const IntPairSet& activeLoops, const int nIterations,
	// 	 std::vector<std::vector<std::pair<std::pair<int, int>, double> > > & odoEdgeRelate2LC_Error)
	IntPairDoubleMap optimize(const IntPairSet& activeLoops, const int nIterations,
		std::vector<std::pair<std::pair<int, int>, double> > & odoEdgeRelateLC_Error,
		std::vector<std::array<double, 5>> & odoEdgeError)
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;
		
		odoEdgeRelateLC_Error.clear();
		odoEdgeError.clear();

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;
			activeEdges.insert( loopclosureEdges[*it]);
		}



		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);
		optimizer->findGauge()->setFixed(true);
		///////////////////
		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
		//optimizer->vertex(0)->setFixed(true);
		optimizer->initializeOptimization(activeEdges);
		optimizer->optimize(nIterations,false);
		optimizer->computeActiveErrors();

		for(int i = 0; i < odoSize; i++)
		{
			// odoEdgeError.push_back(dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
			element_error[0]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
			element_error[1]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(0);
			element_error[2]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(1);
			element_error[3]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(2);
			information_Angle = (dynamic_cast< EdgeType* >(odomteryEdges[i])->information())(2,2);
			element_error[4]  = element_error[3]*element_error[3]*information_Angle;
			odoEdgeError.push_back(element_error);
		}


		double sumLoopChieErr = 0;
		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
			loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();
			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];

			//calculate the error of the odo neighboures of the loop
			//judge the relationship between the node serial and the odometry size

			// node = (*it).first;
			// endNode = (*it).second;	
			// for(int i=0; i<2; i++)
			// {
			// 	if(i == 1)
			// 		node = endNode;
			// 	if(node == 0)
			// 	{
	 	// 			odoErrorElement.first  = std::pair<int,int>(0,0);
	 	// 			odoErrorElement.second = -1;
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);

	  // 				odoErrorElement.first  = std::pair<int,int>(0,1);
	 	// 			odoErrorElement.second = dynamic_cast< EdgeType* >(odomteryEdges[0])->chi2();
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);				
			// 	}
			// 	else if(node == odoSize)
			// 	{
	 	// 			odoErrorElement.first  = std::pair<int,int>(node-1,node);
	 	// 			odoErrorElement.second = dynamic_cast< EdgeType* >(odomteryEdges.back())->chi2();
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);

	  // 				odoErrorElement.first  = std::pair<int,int>(node,node);
	 	// 			odoErrorElement.second = -1;
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);	
			// 	}
			// 	else
			// 	{
	 	// 			odoErrorElement.first  = std::pair<int,int>(node-1,node);
	 	// 			odoErrorElement.second = dynamic_cast< EdgeType* >(odomteryEdges[node-1])->chi2();
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);

	  // 				odoErrorElement.first  = std::pair<int,int>(node,node+1);
	 	// 			odoErrorElement.second = dynamic_cast< EdgeType* >(odomteryEdges[node])->chi2();
	 	// 			odoEdgeRelateLC_Error.push_back(odoErrorElement);	
			// 	}				
			// }
		}

		// NOTE : The number of edges involved is the last element
		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
		//
		return loopClosureLinkError;
	}


	IntPairDoubleMap optimize_active(const IntPairSet& activeLoops, const int nIterations,
		int & startRegion, int & endRegion, std::vector<std::array<double, 5>> & odoEdgeError,
		string ba, string aa)
	{

	// EdgePtrVector 					odomteryEdges;
	// IntPairToEdgePtrMap				loopclosureEdges;
	
	// g2o::OptimizableGraph::EdgeSet 	activeEdges;
	// g2o::OptimizableGraph::VertexSet 	activeVertexes;
		cout<<" "<<endl;
		cout<<"nodes range form "<<startRegion<<" to "<<endRegion<<endl;
		std::pair<std::pair<int, int>, double> odoErrorElement;
		int node, endNode;
		std::array<double, 5> element_error;
		double information_Angle;

		odoEdgeError.clear();

		activeEdges.clear();
		// activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		activeEdges.insert(odomteryEdges.begin() + startRegion,odomteryEdges.begin() + endRegion);

		// odoEdgeSegment.insert(odomteryEdges.begin() + startRegion,odomteryEdges.begin() + endRegion-1);

		// g2o::OptimizableGraph::EdgeSet::iterator iterBegin = odoEdgeSegment.begin(), iterEnd = odoEdgeSegment.end();
		// for(; iterBegin != iterEnd; iterBegin++ )
		// {

		// }

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;
			activeEdges.insert( loopclosureEdges[*it]);
		}

		IntPairDoubleMap loopClosureLinkError;
		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);

		// optimizer->findGauge()->setFixed(true);

		///////////////////
		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
		//optimizer->vertex(0)->setFixed(true);

		optimizer->vertex(endRegion)->setFixed(true);

		// optimizer->vertices().clear();

		cout<<"vertex count before initialization: "<<vertexCount()<<endl;
		cout<<"edgeCount  before initialization: "<<edgeCount()<<endl;
		cout<<"activeEdges.size: "<<optimizer->activeEdges().size()<<endl;
		optimizer->initializeOptimization(activeEdges);
		cout<<"vertex count after initialization: "<<vertexCount()<<endl;
		cout<<"edgeCount  before initialization: "<<edgeCount()<<endl;
		cout<<"activeEdges.size: "<<optimizer->activeEdges().size()<<endl;
		// std::cin.get();


		// const char *g2f=ba.data();
		// std::ofstream out(g2f);
		// optimizer->saveSubset(out,activeEdges);
		// cout<<"active chi2 in saveOptimized Result : "<<optimizer->activeChi2()<<endl; 
		// out.close();

		optimizer->optimize(nIterations,false);

		optimizer->computeActiveErrors();

		optimizer->vertex(endRegion)->setFixed(false);

	// // store();
	// 	const char *g2f___=aa.data();
	// 	std::ofstream out2(g2f___);
	// 	optimizer->saveSubset(out2,activeEdges);
	// 	cout<<"active chi2 in saveOptimized Result : "<<optimizer->activeChi2()<<endl; 
	// 	out2.close();

	// 	cout<<"have sve file after active optimization"<<endl;
	// 	exit(0);

		// for(int i = 0; i < odoSize; i++)
		// {
		// 	// odoEdgeError.push_back(dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
		// 	element_error[0]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->chi2());
		// 	element_error[1]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(0);
		// 	element_error[2]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(1);
		// 	element_error[3]  = (dynamic_cast< EdgeType* >(odomteryEdges[i])->error())(2);
		// 	information_Angle = (dynamic_cast< EdgeType* >(odomteryEdges[i])->information())(2,2);
		// 	element_error[4]  = element_error[3]*element_error[3]*information_Angle;
		// 	odoEdgeError.push_back(element_error);
		// }


		double sumLoopChieErr = 0;
		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			//dynamic_cast< EdgeType* >(loopclosureEdges[*it])->computeError();
			loopClosureLinkError[*it] = dynamic_cast< EdgeType* >(loopclosureEdges[*it])->chi2();
			sumLoopChieErr = sumLoopChieErr + loopClosureLinkError[*it];

		}

		// NOTE : The number of edges involved is the last element
		loopClosureLinkError[IntPair(-1,0)] = optimizer->activeChi2(); // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-1,-1)] = optimizer->activeEdges().size();

		if(loopClosureLinkError[IntPair(-1,0)] > 1)
		{
			cout<<"active optimize chi sum error is large than 1"<<endl;
			// std::cin.get();
		}

		loopClosureLinkError[IntPair(-2,0)] = sumLoopChieErr; // There can be no links with negative IDS
		loopClosureLinkError[IntPair(-2,-1)] = activeLoops.size();		
		//
		return loopClosureLinkError;
	}

	bool write(const char* filename, const IntPairSet& correctLoops)
	{
		restore();
		std::ofstream out(filename);

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for ( IntPairSet::const_iterator it = correctLoops.begin(), end = correctLoops.end();
				it!=end;
				it++)
		{
			activeEdges.insert(loopclosureEdges[*it]);
		}
		optimizer->saveSubset(out,activeEdges);
		out.close();

		return true;
	}
	// optimizer.save("../data/sphere_after.g2o");
	bool saveOptimizedResult(const char* filename, IntPairSet & activeLoops)
	{
		std::ofstream out(filename);

		activeEdges.clear();
		activeEdges.insert(odomteryEdges.begin(),odomteryEdges.end());

		for(
			IntPairSet::const_iterator it = activeLoops.begin(), end = activeLoops.end();
			it!=end;
			it++)
		{
			// if(*it == (std::pair<int,int> (4190, 4085)))
			// 	continue;
			activeEdges.insert( loopclosureEdges[*it]);
		}

		restore();
		//
		optimizer->setVerbose(false);
	    //optimizer->setVerbose(true);
		optimizer->findGauge()->setFixed(true);
		///////////////////
		//std::cerr<<"optimizer->findGauge() output: "<<(*(optimizer->findGauge()))<<std::endl;;
		//optimizer->vertex(0)->setFixed(true);
		optimizer->initializeOptimization(activeEdges);
		optimizer->optimize(5,false);
		optimizer->computeActiveErrors();

		store();
		// optimizer->save(filename);
		optimizer->saveSubset(out,activeEdges);
		cout<<"active chi2 in saveOptimized Result : "<<optimizer->activeChi2()<<endl; 
		out.close();
		return true;
	}

	int edgeDimension()
	{
		return EdgeType::Dimension;
	}

	bool isInitialized()
	{
		return initialized;
	}

	bool removeEdges(const IntPairSet& falseLinks)
	{
		for (IntPairSet::const_iterator it = falseLinks.begin(), end = falseLinks.end(); it!=end; ++it)
		{
			if(optimizer->removeEdge(loopclosureEdges[*it]))
			{
				loopclosureEdges.erase(*it); // Clear up from our records as well
			}
		}
		return true;
	}
};


#endif /* BACKEND_G2O_HPP_ */
