/*
 * test.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: yasir
 */

#include "include/RRR.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <sys/time.h>
#include <assert.h>

#include "g2o/types/slam2d/edge_se2.h"
#include "g2o/types/slam2d/vertex_se2.h"
#include "g2o/types/slam3d/edge_se3.h"
#include "g2o/types/slam3d/vertex_se3.h"
#include "g2o/core/eigen_types.h"

#include <string>


typedef RRR < G2O_Interface
				<
				g2o::VertexSE2, g2o::EdgeSE2
				>
			>
			RRR_2D_G2O;

typedef RRR < G2O_Interface
				<
				g2o::VertexSE3, g2o::EdgeSE3
				>
			>
			RRR_3D_G2O;


/*!
 * This example show how to apply RRR to already existing instance of g2o optimizer.
 * To simulate this, a graph is read into the optimizer and then a pointer of this
 * optimizer is passed to an instance of RRR.
 *
 *
 */
using namespace Eigen;

int main(int argc, char** argv)
{

	if(argc < 2)
	{
		std::cerr<<"Please specify a graph file to read " <<std::endl;
		std::cerr<<argv[0]<<" graph_file clusteringThreshold[default=50]"<<std::endl;
		return -1;
	}

	// std::cin.get();
	// assert(5 == 7);

	double clusteringThreshold = 1;
	int nIter = 4;
	int num_all_vertex = 0;
	string  fullName;

	struct timeval t1, t2, clusterStart, clusterEnd;



	if(argc == 3)
	{
	  clusteringThreshold = atof(argv[2]);
	  // cout<<"clusteringThreshold: "<<clusteringThreshold<<endl;
	}


	std::cout<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<"  Robust SLAM based on Consistent clustering "<<std::endl;
	std::cout<<"   Zihao Xu, Cesar Cadena , Roland Siegwart  "<<std::endl;
	std::cout<<"                ASL, ETH 2018                "<<std::endl;
	std::cout<<"---------------------------------------------"<<std::endl;
	std::cout<<std::endl;
	
	
	std::cout<<"Reading from file :"<<argv[1]<<std::endl;
	if(argc == 3)
	{
		std::cout<<"Clustering threshold (t_g) :" << clusteringThreshold<<std::endl;
	}

	fullName = argv[1];
	int posfile = fullName.find_last_of('/');
	string path = fullName.substr(0, posfile+1); 
	cout<<"path: "<<path<<endl;
	// exit(0);
	string  fileName=fullName.substr(posfile+1);
		cout<<"fileName: "<<fileName<<endl;

    fileName = fileName.substr(0, fileName.length()-4);
	cout<<"fileName"<<fileName <<endl;

	string g2ofile,resultfile;

	if(argc == 3)
	{
		g2ofile = path+"resolved_"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
		resultfile = path+"result_"+fileName+'_'+to_string(int(clusteringThreshold*100))+".g2o";
	}
	else
	{
		g2ofile = path+"resolved_"+fileName+'_'+".g2o";
		resultfile = path+"result_"+fileName+'_'+".g2o";
	}

    cout<<"g2ofile:"<<g2ofile<<endl;
    cout<<"path+fileName: "<<resultfile<<endl;

	
    const char *g2f=g2ofile.data();
    const char *resultFile=resultfile.data();

	// exit(0);
	/* Allocate a g2o optmizer :
	 * It is important that the solver be define in case we want to
	 * reason on a graph already in memory.
	 * */

	g2o::SparseOptimizer optimizer;

	SlamLinearSolver* linearSolver = new SlamLinearSolver();
	linearSolver->setBlockOrdering(false);
	SlamBlockSolver* blockSolver = new SlamBlockSolver(std::unique_ptr<SlamLinearSolver>(linearSolver));
	g2o::OptimizationAlgorithmGaussNewton* solverGauss   = new g2o::OptimizationAlgorithmGaussNewton(std::unique_ptr<SlamBlockSolver>(blockSolver));


	optimizer.setAlgorithm(solverGauss);

	/* load the graph file in the optimizer */
	optimizer.load(argv[1]);


	/* Initialized RRR with the parameters defined */
	RRR_2D_G2O rrr(clusteringThreshold, nIter, fullName);
	// RRR_2D_G2O rrr( nIter, fullName);

	/* Pass the current optimizer pointer to rrr */
// struct timeval t1, t2, clusterStart, clusterEnd;
	gettimeofday(&clusterStart, NULL);

	rrr.setOptimizer(&optimizer,argv[1]);
	// sleep(1);

	gettimeofday(&clusterEnd, NULL);


	/* Find loop closure that are consistent
	 * If the function is passed a bool variable with value true,
	 * it will automatically elimiate all the wrong loops from the
	 * original optimizer. Otherwise, the function removeIncorrectLoops()
	 * can be called to do the same.
	 */

// struct timeval t1, t2, clusterStart, clusterEnd;
	// gettimeofday(&clusterStart, NULL);

	rrr.robustify();

	gettimeofday(&t2, NULL);
	//那么函数f运行所花的时间为
	double deltaT = (t2.tv_sec-clusterStart.tv_sec)  + (t2.tv_usec-clusterStart.tv_usec)/1000000.0;// 微秒
	cout<<" "<<endl;
	cout<<"total time consumed: "<<deltaT<<endl;
	cout<<" "<<endl;

	deltaT = (clusterEnd.tv_sec-clusterStart.tv_sec) + (clusterEnd.tv_usec-clusterStart.tv_usec) / 1000000.0;// 微秒

	cout<<"clustering time consumed: "<<deltaT<<endl;
	/**
	 * If we didn't remove the wrong loop closures earlier, remove them now
	 */
	rrr.removeIncorrectLoops();

	/*
	 * This will write a graph with only the correct loop closures, even when we have not
	 * elimiated incorrect ones using one of the methods above.
	 * */
	std::cout<<std::endl;
	std::cout<<"Output written to rrr-solved.g2o"<<std::endl;

	// rrr.write(g2f);

	// rrr.write_resolved_result(g2f, resultFile);

	// rrr.saveOptimizedResult(OPre);

	/**
	 * Since we have removed the incorrect ones, this file would be the same
	 * as the one above.
	 */
	
	//optimizer.save("g2o-saved.g2o");


	return 0;
}





