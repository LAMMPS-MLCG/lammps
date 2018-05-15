/*
 * mkdssp.h
 *
 *  Created on: May 8, 2018
 *      Author: birva
 */

#ifndef SRC_MKDSSP_H_
#define SRC_MKDSSP_H_
#include <fstream>

#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#if defined USE_COMPRESSION
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#endif
#include <boost/algorithm/string.hpp>
#include"mas.h"
//#include "lammps.h"
//#include "primitives-3d.h"
//#include "align-2d.h"
#include "pointers.h"
#include "structure.h"

namespace po = boost::program_options;
namespace io = boost::iostreams;
namespace ba = boost::algorithm;
using namespace std;
using namespace LAMMPS_NS;
class mkdssp : protected Pointers
{
	int filecount;

public:
	mkdssp(LAMMPS *);
	int init(int, char**                                                                                                          );
	int reinit();


	MProtein* a;
};





#endif /* SRC_MKDSSP_H_ */
