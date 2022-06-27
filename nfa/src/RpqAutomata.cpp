
#include <iostream>
#include <string>
#include <vector>
#include <regex>
#include "../include/parser.hpp"
#include "../include/regular.hpp"
#include "../include/RpqAutomata.hpp"

using std::pair;
using std::string;
using std::unordered_map;

inline uint64_t fwdTrans(Mask *trans, int current, Mask predB)
{
    return *trans[current] & *predB;
}

inline uint64_t bwdTrans(Mask *trans, int current, Mask predB)
{
    return *trans[current & *predB];
}

inline uint64_t next_1SliceInT(const regularData *nfaData, uint64_t current, uint64_t predicate, int type, Mask b)
{
    if (type == FWD)
        return fwdTrans(nfaData->fwdTrans[0], current, b);
    else
        return bwdTrans(nfaData->bwdTrans[0], current, b);
}

uint64_t RpqAutomata::next(uint64_t current, uint64_t predicate, int type)
{
    const auto mapItem = _B.find(predicate);
    if (mapItem == _B.end())
        return 0;

    const auto bPos = mapItem->second;

    if (_nfaData->slices == 1)
        return next_1SliceInT(_nfaData, current, predicate, type, bPos);

    // ELSE handle multiple slices
    uint64_t next = 0;
    if (type == FWD)
        for (auto i = 0; i < _nfaData->slices; i++)
            next |= fwdTrans(_nfaData->fwdTrans[0], current, bPos);
    else
        for (auto i = 0; i < _nfaData->slices; i++)
            next |= bwdTrans(_nfaData->bwdTrans[0], current, bPos);

    return next;
}

regularData *buildNFA(const char *pattern, int length)
{
    Tree *tree;
    Tree **pos;
    int m;

    m = length;
    pos = new Tree *[m];
    tree = parse(pattern, m, pos);

    if ((tree == NULL))
    {
        freeTree(tree);
        free(pos);
        return NULL;
    }

    return regularPreproc(pattern, tree, pos);
}

pair<string, unordered_map<uint64_t, char>> rpqToStringPattern(const string &rpq, const unordered_map<string, uint64_t> &predicates)
{
    unordered_map<uint64_t, char> map;
    char firstUnassigned = 1;
    auto pattern = rpq;

    std::regex predRE("(<[^>]*>)");
    std::sregex_iterator next(rpq.begin(), rpq.end(), predRE);
    std::sregex_iterator end;

    while (next != end)
    {
        std::smatch match = *next;
        const auto str = match.str();

        const auto predicateId = predicates.at(str);
        if (map.find(predicateId) == map.end())
        {
            map[predicateId] = firstUnassigned;

            std::regex replaceRE(str);
            pattern = std::regex_replace(pattern, replaceRE, string(1, firstUnassigned));

            ++firstUnassigned;
        }
        next++;
    }

    std::regex slashRE("/");
    pattern = std::regex_replace(pattern, slashRE, "");

    return {pattern, map};
}

pair<regularData *, unordered_map<uint64_t, char>> build(const string &rpq, const unordered_map<string, uint64_t> &predicates)
{

    // 1. Transform RPQ to a regular expression
    auto patternAndMap = rpqToStringPattern(rpq, predicates);
    auto pattern = patternAndMap.first;

    // 2. Use nrgrep to parse the pattern and build the NFA
    auto regularData = buildNFA(pattern.c_str(), pattern.length());

    return {regularData, patternAndMap.second};
}

pair<regularData *, unordered_map<uint64_t, char>> build(const string &pattern, const unordered_map<uint64_t , char> &map)
{


    // 2. Use nrgrep to parse the pattern and build the NFA
    auto regularData = buildNFA(pattern.c_str(), pattern.length());

    return {regularData, map};
}

RpqAutomata::RpqAutomata() {
    _nfaData = NULL;
    _isValid = false;
};

RpqAutomata::RpqAutomata(const string &rpq, const unordered_map<string, uint64_t> &predicates)
{
    auto nfaAndMap = build(rpq, predicates);

    _nfaData = nfaAndMap.first;
    _isValid = _nfaData != NULL;

    for (auto pair : nfaAndMap.second)
    {
        _B[pair.first] = _nfaData->B[pair.second];
    }
}

RpqAutomata::RpqAutomata(const string &pattern, const unordered_map<uint64_t , char> &map)
{
    auto nfaAndMap = build(pattern, map);

    _nfaData = nfaAndMap.first;
    _isValid = _nfaData != NULL;

    for (auto pair : nfaAndMap.second)
    {
        _B[pair.first] = _nfaData->B[pair.second];
    }
}

unordered_map<uint64_t, uint64_t> RpqAutomata::getB()
{
    unordered_map<uint64_t, uint64_t> result;

    for (auto pair : _B)
    {
        result[pair.first] = *pair.second;
    }

    return result;
}

RpqAutomata::~RpqAutomata()
{
    if (_nfaData != NULL)
    {
        regularFree(_nfaData);
        _nfaData = NULL;
    }
}