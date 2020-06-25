#include "node.h"

Node::Node():
effMass( 0.0 )
{}
Node::~Node()
{

	memPtr->Delete2ndMat( coord, npoin );
	memPtr->Delete2ndMat( cordi, npoin );
	memPtr->Delete2ndMat( coordUnWrapped, npoin );
	memPtr->Delete2ndMat( cordiUnWrapped, npoin );
	delete [] dispt;
	delete [] dispt0;
	delete [] velot;
	delete [] veloi;
	delete [] fintl;
	delete [] absfl;
	delete [] mass;
	delete [] dldis;
	delete [] presVeloc;
	delete [] fixID;
	delete [] rigidID;
	delete [] extForce;
}
