#include "BeliefTreeNode.h"

BeliefTreeNode::BeliefTreeNode(void)
{
	s = NULL;
}

BeliefTreeNode::~BeliefTreeNode(void)
{
}
void BeliefTreeNode::print()
{
	cout << "BeliefTreeNode:" << endl;
	cout << "CacheIndex row sval " << cacheIndex.row << " " << cacheIndex.sval << endl;
	cout << "s belief ";
	s->bvec->write(cout) << endl;    // calls SparseVector::write(std::ostream& out)
	cout << " sval " << s->sval << endl;
	//cout << " addresses " << &(s->sval) << " addresses " << &(cacheIndex.sval) << endl << endl;

	for(vector<BeliefTreeQEntry>::iterator iter1 = Q.begin() ; iter1 != Q.end() ; iter1 ++)
	{
		BeliefTreeQEntry&  entry= *iter1;
		for(vector<BeliefTreeObsState*>::iterator iter2 = entry.stateOutcomes.begin() ; iter2 != entry.stateOutcomes.end() ; iter2 ++)
		{
			BeliefTreeObsState* obsState = *iter2;
			if( obsState == NULL)
			{
				continue;
			}

			for(vector<BeliefTreeEdge*>::iterator iter3 = obsState->outcomes.begin() ; iter3 != obsState->outcomes.end() ; iter3 ++)
			{
				BeliefTreeEdge* edge = *iter3;
				if(edge!= NULL)
				{
					edge->nextState->print();
				}
			}

		}
	}
}



BeliefTreeEdge::BeliefTreeEdge()
{
	nextState = NULL;
}
BeliefTreeEdge::~BeliefTreeEdge()
{

}
