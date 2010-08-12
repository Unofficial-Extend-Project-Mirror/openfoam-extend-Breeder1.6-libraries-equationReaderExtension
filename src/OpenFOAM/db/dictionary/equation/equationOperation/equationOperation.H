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

Class
    Foam::equationOperation

Description
    Defines a single operation to be performed in sequence while evaluating an
    equation read from a dictionary.

SourceFiles
    equationOperationI.H
    equationOperation.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationOperation_H
#define equationOperation_H

#include "word.H"
#include "label.H"
#include "IOstreams.H"
// #include "Istream.H"
// #include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class equationOperation;

// Friend Operators

int operator==(const equationOperation&, const equationOperation&);
int operator!=(const equationOperation&, const equationOperation&);

// Friend IOstream Operators

Istream& operator>>(Istream&, equationOperation&);
Ostream& operator<<(Ostream&, const equationOperation&);



/*---------------------------------------------------------------------------*\
                        Class equationOperation Declaration
\*---------------------------------------------------------------------------*/

class equationOperation
{

public:

        enum sourceListType
        {
            slnone,
            sldictSource,
            slexternalDScalar,
            slexternalScalar,
            slexternalScalarList,
            slinternalScalar,
            slequation,
            slstorage
        };

        enum operationType
        {
            otnone,
            otretrieve,
            otstore,
            otplus,
            otminus,
            ottimes,
            otdivide,
            otpow,
            otsign,
            otpos,
            otneg,
            otmag,
            otlimit,
            otminMod,
            otsqrtSumSqr,
            otsqr,
            otpow3,
            otpow4,
            otpow5,
            otinv,
            otsqrt,
            otcbrt,
            othypot,
            otexp,
            otlog,
            otlog10,
            otsin,
            otcos,
            ottan,
            otasin,
            otacos,
            otatan,
            otsinh,
            otcosh,
            ottanh,
            otasinh,
            otacosh,
            otatanh,
            oterf,
            oterfc,
            otlgamma,
            otj0,
            otj1,
            otjn,
            oty0,
            oty1,
            otyn,
            otmax,
            otmin,
            otstabilise
        };

private:

    // List to read the data from
    sourceListType sourceList_;
    
    // Index in the list where the data is located.  Note, the
    // equationOperationLists created by equationReader make this variable a
    // 1-indexed (i.e. starts from 1, not zero) in order to use its sign to
    // store the sign of the variable.  The sourceLists are zero-indexed, so
    // equationReader will constantly be adding / subtracting 1 to get these
    // to match
    label sourceIndex_;
    
    // Rather than store the keywords that have to be searched in a dictionary,
    // equationReader keeps its own list of keywords, and the dictLookupIndex
    // is the index in this list.  This is only applicable if the sourceList is
    // of type sldictSource
    label dictLookupIndex_;
    
    // The operation to be performed (+ - sin exp min, etc...)
    operationType operation_;

public:

    // Static data members
    
        static const char* const typeName;

    // Constructors

        //- Construct null
        equationOperation();
        
        //- Construct from components
        equationOperation
        (
            sourceListType sourceList,
            label sourceIndex,
            label dictLookupIndex,
            operationType operation
        );


    // Destructor
    ~equationOperation();

    // Member functions

        // Access
        
            //- Const access to source list
            inline const sourceListType& sourceList() const;
            
            //- Access to source list
            inline sourceListType& sourceList();

            //- Const access to source index
            inline const label& sourceIndex() const;
            
            //- Access to source index
            inline label& sourceIndex();

            //- Const access to dictionary lookup name index
            inline const label& dictLookupIndex() const;
            
            //- Access to dictionary lookup name index
            inline label& dictLookupIndex();

            //- Const access to operation
            inline const operationType& operation() const;
            
            //- Access to operation
            inline operationType& operation();
            
            //- Look up operation number, given a word
            static operationType findOp(const word& opName);

            //- Look up operation name
            static word opName(const operationType& op);

            //- Look up sourceList name
            static word sourceName(const sourceListType& sl);

    
    // Friend IOstream Operators

        friend Istream& operator>>(Istream&, equationOperation&);
        friend Ostream& operator<<(Ostream&, const equationOperation&);
      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOperationI.H"

#endif
