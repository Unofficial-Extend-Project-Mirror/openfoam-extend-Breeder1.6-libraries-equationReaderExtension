/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "equationReader.H"
#include "mathematicalConstants.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//const char* const Foam::equationReader::typeName = "equationReader";
    defineTypeNameAndDebug(equationReader, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::equationReader::fatalParseError
(
    const Foam::label index,
    const Foam::tokenList& tl,
    const Foam::label fromToken,
    const Foam::label toToken,
    const string& errorIn,
    const OStringStream& description
)
{
    OStringStream errorMessage;
    forAll(tl, i)
    {
        if (i == fromToken)
        {
            errorMessage << "<";
        }
        if (tl[i].isPunctuation() && tl[i].pToken() == token::COLON)
        {
            errorMessage << "^";
        }
        else
        {
            errorMessage << tl[i];
        }
        if (i == toToken)
        {
            errorMessage << ">";
        }
    }

    FatalErrorIn(errorIn) << "Parsing error in the equation for "
        << eqns_[index].equationName() << ", given by:" << endl
        << endl << token::TAB << eqns_[index].rawText() << endl
        << endl << "Error occurs withing the < angle brackets >:" << endl
        << endl << token::TAB << errorMessage.str() << endl << endl
        << description.str()
        << abort(FatalError);
}


Foam::string Foam::equationReader::stringPreconditioner(const string& rawText)
{
    string rawWorking(rawText);

    // Negative exponent workaround
    for (label i = 0; i < 10; i++)
    {
        string strTemp(name(i));
        string strTemp2(name(i));
        strTemp.append("e-");
        strTemp2.append("&");
        stringReplaceAll(rawWorking, strTemp, strTemp2);
    }
    
    stringReplaceAll(rawWorking, "^", " : ");
    stringReplaceAll(rawWorking, "(", " ( ");
    stringReplaceAll(rawWorking, ")", " ) ");
    stringReplaceAll(rawWorking, "+", " + ");
    stringReplaceAll(rawWorking, "-", " - ");
    stringReplaceAll(rawWorking, "&", "e-");
    stringReplaceAll(rawWorking, "*", " * ");
    stringReplaceAll(rawWorking, "/", " / ");
    stringReplaceAll(rawWorking, ",", " , ");
    return rawWorking;
}


void Foam::equationReader::stringReplaceAll
(
    Foam::string& working,
    const Foam::string& findMe,
    const Foam::string& replaceWith
)
{
    size_t position(0);
    size_t offset(0);
    string subWorking(working);
    while (position != string::npos)
    {
        position = subWorking.string::find(findMe);
        if (position != string::npos)
        {
            working.::std::string::replace
            (
                position + offset, findMe.size(), replaceWith
            );
            offset += position + replaceWith.size();
            subWorking = subWorking.substr(position + findMe.size());
        }
    }
}


Foam::labelList Foam::equationReader::findMaxParenthesis
(
    const Foam::labelList& parenthesisList,
    const Foam::labelList& equationIndices
) const
{
    labelList returnMe(equationIndices.size());
    label currentMax(-1);
    label atIndex(0);
    bool groupDone(false);

    forAll(equationIndices, i)
    {
        if (mag(parenthesisList[equationIndices[i]]) > currentMax)
        {
            groupDone = false;
            atIndex = 0;
            currentMax = mag(parenthesisList[equationIndices[i]]);
            returnMe[atIndex] = equationIndices[i];
        }
        else if
        (
            (mag(parenthesisList[equationIndices[i]]) == currentMax)
         && (!groupDone)
        )
        {
            atIndex++;
            returnMe[atIndex] = equationIndices[i];
        }
        else if (mag(parenthesisList[equationIndices[i]]) < currentMax)
        {
            groupDone = true;
        }
    }
    returnMe.setSize(atIndex + 1);
    return returnMe;
}


Foam::labelList Foam::equationReader::findMaxOperation
(
    const Foam::labelList& opLvl,
    const Foam::labelList& indices
)
{
    label maxOpLvl(-1);
    label atIndex(-1);
    labelList returnMe(indices.size());
    bool groupDone(false);

    forAll(indices, i)
    {
        if (opLvl[indices[i]] > maxOpLvl)
        {
            atIndex = 0;
            returnMe[0] = indices[i];
            if (i > 0)
            {
                atIndex++;
                returnMe[0] = indices[i - 1];
                returnMe[1] = indices[i];
            }
            groupDone = false;
            maxOpLvl = opLvl[indices[i]];
        }
        else if
        (
            (
                (opLvl[indices[i]] == maxOpLvl)
             || (opLvl[indices[i]] == 0)
            )
         && !groupDone
        )
        {
            atIndex++;
            returnMe[atIndex] = indices[i];
        }
        else if ((opLvl[indices[i]] < maxOpLvl) && (opLvl[indices[i]] != 0))
        {
            groupDone = true;
        }
    }
    returnMe.setSize(atIndex + 1);
    return returnMe;
}


void Foam::equationReader::absorbNegatives
(
    const Foam::label equationIndex,
    const Foam::tokenList& tl,
    Foam::labelList& eqnIndices,
    Foam::labelList& subEqnIndices,
    Foam::equationOperationList& map,
    Foam::labelList& opLvl
)
{
    // Negatives are identified by a map with a negative sourceIndex
    // To accommodate this behaviour, source indices are 1-index.
    forAll(subEqnIndices, i)
    {
        if (map[subEqnIndices[i]].dictLookupIndex() == -1)
        {
            if
            (
                (subEqnIndices.size() == i)
             || (opLvl[subEqnIndices[i + 1]] != 0)
            )
            {
                OStringStream description;
                description << "Misplaced negative / subtraction operator.";
                fatalParseError
                (
                    equationIndex,
                    tl,
                    subEqnIndices[i],
                    subEqnIndices[i],
                    "equationReader::absorbNegatives",
                    description
                );
            }
            map[subEqnIndices[i + 1]].sourceIndex() = 
                -map[subEqnIndices[i + 1]].sourceIndex();
            
            trimListWithParent(eqnIndices, subEqnIndices, i, i);
        }
    }
}


void Foam::equationReader::dsEqual
(
    Foam::dimensionedScalar& dso,
    const Foam::dimensionedScalar& dsi
)
{
    dso.value() = dsi.value();
    dso.name() = dsi.name();
    dso.dimensions().reset(dsi.dimensions());
}


void Foam::equationReader::trimListWithParent
(
    Foam::labelList& parent,
    Foam::labelList& indices,
    Foam::label from,
    Foam::label to,
    Foam::label exceptFor
)
{
    if
    (
        (to > (indices.size() - 1))
     || (from < 0)
     || (from > (indices.size() - 1))
     || (to < from)
    )
    {
        FatalErrorIn
        (
            "equationReader::trimListWithParent(parent, indices, from, to, "
            "exceptFor)"
        )
            << "Bad indices.  parent is " << parent << ", indices are "
            << indices << ", from is " << from << ", to is " << to
            << " exceptFor is " << exceptFor << "."
            << abort(FatalError);
    }

    for (label i(from); i <= to; i++)
    {
        label removeMe(indices[i]);
        if (i == exceptFor) continue;
        forAll(parent, j)
        {
            if (parent[j] == removeMe)
            {
                trimList(parent, j, j);
                break;
            }
        }
    }
    trimList(indices, from, to, exceptFor);
}


void Foam::equationReader::trimList
(
    Foam::labelList& indices,
    Foam::label from,
    Foam::label to,
    Foam::label exceptFor
)
{
    if
    (
        (to > (indices.size() - 1))
     || (from < 0)
     || (from > (indices.size() - 1))
     || (to < from)
    )
    {
        FatalErrorIn
        (
            "equationReader::trimList(indices, from, to, exceptFor"
            "exceptFor)"
        )
            << "Bad indices.  indices are " << indices << ", from is "
            << from << ", to is " << to << " exceptFor is " << exceptFor << "."
            << abort(FatalError);
    }

    if (!(exceptFor == from && from == to))
    {
        
        if (exceptFor == from)
        {
            from++;
        }
        else if (exceptFor == to)
        {
            to--;
        }
        else if ((exceptFor < to) && (exceptFor > from))
        {
            indices[from] = indices[exceptFor];
            from++;
        }

        for (label i(from); i < (indices.size() + from - to - 1); i++)
        {
            indices[i] = indices[i + to - from + 1];
        }
        indices.setSize(indices.size() - to + from - 1);
    }
}


Foam::label Foam::equationReader::findIndex
(
    const Foam::label value,
    const Foam::labelList& indexList
) const
{
    forAll (indexList, i)
    {
        if (indexList[i] == value)
        {
            return i;
        }
    }
    return -1;
}


Foam::equationOperation Foam::equationReader::findSource
(
    const Foam::word& varName,
    const Foam::label equationIndex
)
{
    // Search order:
    // -other equations
    // -externalDScalars
    // -externalScalars
    // -dictSources
    equationOperation returnMe
    (
        equationOperation::slnone,
        0,
        0,
        equationOperation::otnone
    );

    // Searching known equations
    for (label eqs(0); eqs < eqns_.size(); eqs++)
    {
        if (eqns_[eqs].equationName() == varName)
        {
            returnMe = equationOperation
            (
                equationOperation::slequation,
                eqs + 1,
                0,
                equationOperation::otnone
            );
        }
    }
    
    forAll(externalDScalars_, i)
    {
        if (externalDScalars_[i].name() == varName)
        {
            returnMe = equationOperation
            (
                equationOperation::slexternalDScalar,
                i + 1,
                0,
                equationOperation::otnone
            );
        }
    }
    
    if (externalScalars_.size() != externalScalarNames_.size())
    {
        FatalErrorIn("equationReader::findSource")
            << "Size mismatch detected in externalScalars = "
            << externalScalars_.size() << " and externalScalarNames = "
            << externalScalarNames_.size() << "."
            << abort(FatalError);
    }
    
    forAll(externalScalars_, j)
    {
        if (externalScalarNames_[j] == varName)
        {
            returnMe = equationOperation
            (
                equationOperation::slexternalScalar,
                j + 1,
                0,
                equationOperation::otnone
            );
        }
    }

    forAll(dictSources_, i)
    {
        if (dictSources_[i].found(varName, true))
        {
            label dictLookupIndex(-1);
//            wordList& wl(eqns_[equationIndex].dictLookups());

            ITstream srcStrm
            (
                dictSources_[i].lookup(varName, true)
            );
            if (isDimensionedScalar(srcStrm) || isScalar(srcStrm))
            {
                // Look for varName in dictLoookup names, append if not found
                forAll(dictLookups_, j)
                {
                    if (varName == dictLookups_[j])
                    {
                        dictLookupIndex = j;
                        break;
                    }
                }
                if (dictLookupIndex < 0)
                {
                    dictLookups_.setSize(dictLookups_.size() + 1);
                    dictLookups_.set
                    (
                        dictLookups_.size() - 1,
                        new word(varName)
                    );
                    dictLookupIndex = dictLookups_.size() - 1;
                }
                
                returnMe = equationOperation
                (
                    equationOperation::sldictSource,
                    i + 1,
                    dictLookupIndex,
                    equationOperation::otnone
                );
            }
            else if (isEquation(srcStrm))
            {
                // Is it a known equation already?
                for (label eqs(0); eqs < eqns_.size(); eqs++)
                {
                    if (eqns_[eqs].equationName() == varName)
                    {
                        returnMe = equationOperation
                        (
                            equationOperation::slequation,
                            eqs + 1,
                            0,
                            equationOperation::otnone
                        );
                        break;
                    }
                }
                
                if (returnMe.sourceList() == equationOperation::slnone)
                {
                    // The equation has not been read yet.  Create an unparsed
                    // equation.  It will be parsed during evaluate, or if it
                    // is later read.
                    equation eqn(srcStrm);
                    eqn.equationName() = varName;
                    createEquation(eqn);
/*
                    token it(srcStrm);
                    createEquation
                    (
                        equation
                        (
                            varName,
                            it.stringToken()
                        )
                    );
*/
                    returnMe = equationOperation
                    (
                        equationOperation::slequation,
                        eqns_.size(),
                        0,
                        equationOperation::otnone
                    );
                    break;
                }
            }
            break;
        } // end if varName is in dictionary
    } // end dictionary search loop
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSource
(
    const Foam::label equationIndex,
    const Foam::label equationOperationIndex,
    const Foam::label maxStoreIndex,
    const Foam::label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    switch (eqOp.sourceList())
    {
        case equationOperation::slnone:
            dsEqual(returnMe, dimensionedScalar("noSource", dimless, 0));
            break;
        case equationOperation::sldictSource:
        {
            if (zeroSourceIndex >= dictSources_.size())
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << dictSources_.size() - 1 << ")"
                    << abort(FatalError);
            }

            // Accepts scalars, dimensionedScalars
            word varName(dictLookups_[eqOp.dictLookupIndex()]);
            
            ITstream srcStrm
            (
                dictSources_[zeroSourceIndex].lookup(varName, true)
            );
            if (isDimensionedScalar(srcStrm))
            {
                srcStrm >> returnMe;
            }
            else if (isScalar(srcStrm))
            {
                returnMe.name() = varName;
                returnMe.value() = readScalar(srcStrm);
            }
            else
            {
                // Neither scalar nor dimensionedScalar
                FatalIOErrorIn
                (
                    "equationReader::getSource",
                    dictSources_[zeroSourceIndex]
                )
                    << "Expecting a scalar or a dimensionedScalar.  Keyword "
                    << varName << " is referenced by an equation, and therfore"
                    << " can only be one of these two types."
                    << exit(FatalIOError);
            }
            break;
        }
        case equationOperation::slexternalDScalar:
            if (zeroSourceIndex >= externalDScalars_.size())
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << externalDScalars_.size() - 1 << ")"
                    << abort(FatalError);
            }
            dsEqual(returnMe, externalDScalars_[zeroSourceIndex]);
            break;
        case equationOperation::slexternalScalar:
            if (zeroSourceIndex >= externalScalars_.size())
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << externalScalars_.size() - 1 << ")"
                    << abort(FatalError);
            }
            returnMe.name() = externalScalarNames_[zeroSourceIndex];
            returnMe.value() = externalScalars_[zeroSourceIndex];
            break;
        case equationOperation::slinternalScalar:
            if (zeroSourceIndex >= internalScalars_.size())
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << internalScalars_.size() - 1 << ")"
                    << abort(FatalError);
            }
            returnMe.name() = "internalConstant";
            returnMe.value() =
                internalScalars_[zeroSourceIndex];
            break;
        case equationOperation::slequation:
            if  (zeroSourceIndex >= eqns_.size())
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << eqns_.size() - 1 << ")"
                    << abort(FatalError);
            }
            
            // Check for circular references
            dependents_.setSize(dependents_.size() + 1);
            dependents_[dependents_.size() - 1] = equationIndex;
            forAll(dependents_, i)
            {
                if (dependents_[i] == zeroSourceIndex)
                {
                    // Circular reference detected
                    
                    string dependencies;
                    for (label j(i); j < dependents_.size(); j++)
                    {
                        dependencies.append
                        (
                            eqns_[j].equationName()
                        );
                        dependencies.append("-->");
                    }
                    dependencies.append(eqns_[i].equationName());
                    FatalErrorIn("masterEquation::getSource")
                        << "Circular reference detected when evaluating "
                        << "the equation for " << eqn.equationName()
                        << ", given by:" << token::NL << token::TAB
                        << eqn.rawText() << token::NL << "The circular "
                        << "dependency is:" << token::NL << token::TAB
                        << dependencies
                        << abort(FatalError);
                }
            }
            if (debug)
            {
                Info << "Embedded equation dispatch." << endl;
            }
            dsEqual(returnMe, evaluate(zeroSourceIndex, maxStoreIndex + 1));
            if (debug)
            {
                Info << "Returned from embedded equation." << endl;
            }
            break;
        case equationOperation::slstorage:
            if ((zeroSourceIndex + storageOffset) > maxStoreIndex)
            {
                FatalErrorIn("equationReader::getSouce")
                    << "Index " << zeroSourceIndex << " out of bounds (0, "
                    << maxStoreIndex - storageOffset << ")"
                    << abort(FatalError);
            }
            dsEqual(returnMe, storage_[zeroSourceIndex + storageOffset]);
            break;
    }
    
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::label Foam::equationReader::addInternalScalar(const Foam::scalar& value)
{
    forAll(internalScalars_, i)
    {
        if (mag(internalScalars_[i] - value) < VSMALL)
        {
            return i;
        }
    }
    internalScalars_.setSize(internalScalars_.size() + 1);
    internalScalars_.set(internalScalars_.size() - 1, new scalar(value));
    return internalScalars_.size() - 1;
}


bool Foam::equationReader::isDimensionedScalar(Foam::ITstream& is)
{
    tokenList tl(12);
    label found(0);
    while (!is.eof())
    {
        tl[found] = token(is);
        found++;
        if (found > 12)
        {
            is.rewind();
            return false;
        }
    }
    if
    (
        (
            (found == 11)
         && tl[0].isWord()
         && tl[1].isPunctuation()
         && tl[2].isNumber()
         && tl[3].isNumber()
         && tl[4].isNumber()
         && tl[5].isNumber()
         && tl[6].isNumber()
         && tl[7].isNumber()
         && tl[8].isNumber()
         && tl[9].isPunctuation()
         && tl[10].isNumber()
        )
     || (
            (found == 9)
         && tl[0].isWord()
         && tl[1].isPunctuation()
         && tl[2].isNumber()
         && tl[3].isNumber()
         && tl[4].isNumber()
         && tl[5].isNumber()
         && tl[6].isNumber()
         && tl[7].isPunctuation()
         && tl[8].isNumber()
        )
    )
    {
        is.rewind();
        return true;
    }
    else
    {
        is.rewind();
        return false;
    }
}


bool Foam::equationReader::isScalar(ITstream& is)
{
    token firstToken(is);
    if (firstToken.isNumber() && is.eof())
    {
        is.putBack(firstToken);
        return true;
    }
    else
    {
        is.putBack(firstToken);
        return false;
    }
}


bool Foam::equationReader::isEquation(ITstream& is)
{
    token firstToken(is);
    if (firstToken.isString() && is.eof())
    {
        is.putBack(firstToken);
        return true;
    }
    else
    {
        is.putBack(firstToken);
        tokenList tl(12);
        label found(0);
        while (!is.eof())
        {
            tl[found] = token(is);
            found++;
            if (found == 12)
            {
                is.rewind();
                return false;
            }
        }
        if
        (
            (
                (found == 11)
             && tl[0].isWord()
             && tl[1].isPunctuation()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isNumber()
             && tl[8].isNumber()
             && tl[9].isPunctuation()
             && tl[10].isString()
            )
         || (
                (found == 10)
             && tl[0].isPunctuation()
             && tl[1].isNumber()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isNumber()
             && tl[8].isPunctuation()
             && tl[9].isString()
            )
         || (
                (found == 9)
             && tl[0].isWord()
             && tl[1].isPunctuation()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isPunctuation()
             && tl[8].isString()
            )
         || (
                (found == 8)
             && tl[0].isPunctuation()
             && tl[1].isNumber()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isPunctuation()
             && tl[7].isString()
            )
        )
        {
            is.rewind();
            return true;
        }
        else
        {
            is.rewind();
            return false;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationReader::equationReader()
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationReader::~equationReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationReader::addDataSource(const Foam::dictionary * dict)
{
    dictSources_.setSize(dictSources_.size() + 1);
    dictSources_.set(dictSources_.size() - 1, dict);
}


void Foam::equationReader::addDataSource
(
    const Foam::scalar * value,
    const Foam::word name
)
{
    externalScalars_.setSize(externalScalars_.size() + 1);
    externalScalars_.set(externalScalars_.size() - 1, value);

    externalScalarNames_.setSize(externalScalarNames_.size() + 1);
    externalScalarNames_[externalScalarNames_.size() - 1] = name;
}


void Foam::equationReader::addDataSource(const Foam::dimensionedScalar * ds)
{
    externalDScalars_.setSize(externalDScalars_.size() + 1);
    externalDScalars_.set(externalDScalars_.size() - 1, ds);
}


void Foam::equationReader::addDataSource
(
    const Foam::scalarList * scalars,
    const Foam::wordList& names
)
{
    if (scalars->size() != names.size())
    {
        FatalErrorIn("equationReader::addDataSource(const scalarList "
               "scalars, const wordList names)")
            << "Size mismatch.  scalars.size = " << scalars->size()
            << ", names.size() = " << names.size()
            << abort(FatalError);
    }

    for (label i(0); i < scalars->size(); i++)
    {
        addDataSource(&scalars->operator[](i), names[i]);
    }
}


void Foam::equationReader::addDataSource
(
    const Foam::PtrList<dimensionedScalar> * ds
)
{
    for (label i(0); i < ds->size(); i++)
    {
        addDataSource(&ds->operator[](i));
    }
}


void Foam::equationReader::createEquation
(
    Foam::equation eqn,
    Foam::dimensionedScalar * outputVar
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        if (debug)
        {
            Info << "Creating equation " << eqn.equationName() << " at index "
                << eqns_.size() << endl;
        }
        eqns_.setSize(eqns_.size() + 1);
        eqns_.set(eqns_.size() - 1, new equation(eqn));
        outputDScalars_.setSize(eqns_.size());
        outputDScalars_.set(eqns_.size() - 1, outputVar);
    }
    else
    {
        FatalErrorIn("equationReader::createEquation")
            << "Equation already exists.  Use readEquation to suppress error."
            << abort(FatalError);
    }
}


void Foam::equationReader::readEquation
(
    const dictionary * sourceDict,
    const word& eqnName,
    dimensionedScalar * outputVar
)
{
    equation eqn(sourceDict->lookup(eqnName));
    eqn.equationName() = eqnName;

    readEquation
    (
        eqn,
        outputVar
    );
}


void Foam::equationReader::readEquation
(
    Foam::equation eqn,
    Foam::dimensionedScalar * outputVar
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        createEquation(eqn, outputVar);
        parse(eqns_.size() - 1);
    }
    else if
    (
        (eqn.rawText() != eqns_[index].rawText())
     || (!eqns_[index].size())
    )
    {
        parse(index);
        if (outputVar != NULL)
        {
            outputDScalars_.set(index, outputVar);
        }
    }
}


void Foam::equationReader::update(const Foam::word& equationName)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::update(const word)")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    if (outputDScalars_.set(index))
    {
        outputDScalars_[index] = evaluate(index);
    }
}


void Foam::equationReader::update(const Foam::label& index)
{
    if (outputDScalars_.set(index))
    {
        outputDScalars_[index] = evaluate(index);
    }
}


void Foam::equationReader::update()
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        update(i);
    }
}


bool Foam::equationReader::found(const Foam::word& equationName)
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        if (eqns_[i].equationName() == equationName)
        {
            return true;
        }
    }
    return false;
}


Foam::label Foam::equationReader::lookup(const Foam::word& equationName)
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        if (eqns_[i].equationName() == equationName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::equationReader::deleteEquation(const Foam::word& equationName)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        WarningIn("equationReader::deleteEquation(equationName)")
            << "Equation name " << equationName << " not found." << endl;
    }
    deleteEquation(index);
}


void Foam::equationReader::deleteEquation(const Foam::label& index)
{
    if ((index < 0) || (index >= eqns_.size()))
    {
        FatalErrorIn("equationReader::deleteEquation(index)")
            << "Index " << index << " out of bounds (0, " << eqns_.size() - 1
            << ")"
            << abort(FatalError);
    }
    for (label i = index; i < (eqns_.size() - 1); i++)
    {
        eqns_[i] = eqns_[i + 1];
    }
    
    eqns_.setSize(eqns_.size() - 1);
}


void Foam::equationReader::status()
{
    Info << "*** equationReader object status report ***" << endl;
    Info << "Data sources:" << endl;
    Info << token::TAB << dictSources_.size() << " dictionaries, with recorded "
        << "keywords:" << endl;
    forAll(dictLookups_, i)
    {
        Info << token::TAB << token::TAB << dictLookups_[i] << endl;
    }
    Info << token::TAB << externalDScalars_.size() << " external "
        << "dimensionedScalars with names:" << endl;
    forAll(externalDScalars_, i)
    {
        Info << token::TAB << token::TAB << externalDScalars_[i].name()
            << endl;
    }
    Info << token::TAB << externalScalars_.size() << " external "
        << "scalars with names:" << endl;
    forAll(externalScalarNames_, i)
    {
        Info << token::TAB << token::TAB << externalScalarNames_[i] << endl;
    }
    Info << token::TAB << internalScalars_.size() << " internal scalars."
        << endl;
    Info << eqns_.size() << " equations in memory:" << endl;
    forAll(eqns_, i)
    {
        Info << i << ":" << token::TAB << eqns_[i].size() << " operations,"
            << token::TAB << eqns_[i].equationName() << token::TAB;
        if ((outputDScalars_.size() > i) && outputDScalars_.set(i))
        {
            Info << "active output to " << outputDScalars_[i].name() << "."
                << endl;
        }
        else
        {
            Info << "passive output." << endl;
        }
        if (eqns_[i].changeDimensions())
        {
            Info << token::TAB << "Dimension override to: "
                << eqns_[i].overrideDimensions() << endl;
        }
        Info << token::TAB << eqns_[i].rawText() << endl;
    }
}

/*
template <>
dimensioned<scalar>::dimensioned
(
    Istream& is
)
:
    name_(is),
    dimensions_(dimless),
    value_(0)
{
    is.rewind();
    token t(is);

    if (t.isString())
    {
        equationReader eqn;
        eqn.readEquation
        (
            equation
            (
                name_,
                t.stringToken()
            )
        );
        *this = eqn.evaluate(0);
    }
    else
    {
        dimensions_ = dimensionSet(is);
        value_ = scalar(pTraits<scalar>(is));
    }
}
*/

#include "equationReaderCreateMap.C"
#include "equationReaderEvaluate.C"
#include "equationReaderParse.C"

// ************************************************************************* //
