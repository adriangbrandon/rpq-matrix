#ifndef RPQ_AUTOMATA_INCLUDED
#define RPQ_AUTOMATA_INCLUDED

#include <string>
#include <vector>
#include <unordered_map>
#include "regular.hpp"

/* A class abastracting an NFA for an RPQ.
   
   Remarks:
        - Assumes that there are at most 64 different states in the resulting automata
        - Assumes that the predicate names have already been maped to integers, and it 
          requires the map as an argument in the constructor.
        - Predicate names in the RPQ must be enclosed in < > (e.g. <P2119>?/<P2515>).
*/
class RpqAutomata
{

public:
  /* Creates an instance of an RpqAutomata given the RPQ query, and the global map of predicate names to predicate ids
    */
  RpqAutomata(const std::string rpq, const std::unordered_map<std::string, uint64_t> predicates);

   ~RpqAutomata(void);

  /* Returns wether the NFA was correctly constructed from the RPQ
    */
  inline bool isValid() { return _isValid; };

  /* Returns a map with the predicates from the RPQ as keys, and the corresponding row in B as values
    */
  std::unordered_map<uint64_t, uint64_t> getB();

  /* Returns the set of active states after a transition from a given set of states by a given predicate 
       Params:
            - current: a mask representing the set of active states
            - predicate: the predicate used for transitioning to the next states
            - type: can be either FWD or BWD (see basics.h), and indicates whether to
              transition in a backward or forward fashion.
        Returns:
            A mask representing the set of states that are active after the transition
    */
  uint64_t next(uint64_t current, uint64_t predicate, int type);

  /* Returns whether a given set of states contain at least one final state of the automata
    */
  inline bool atFinal(uint64_t activeStates, int type)
  {
    return type == FWD
               ? activeStates & *(_nfaData->final)
               : activeStates & 1;
  }

  inline uint64_t getFinalStates()
  {
    return *(_nfaData->final);
  }

private:
  bool _isValid;
  regularData *_nfaData;
  std::unordered_map<uint64_t, Mask> _B;
};

#endif