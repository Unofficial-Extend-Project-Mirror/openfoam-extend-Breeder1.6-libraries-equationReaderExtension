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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const equationList& equationReader::eqns() const
{
    return eqns_;
}

inline const UPtrList<const dictionary>& equationReader::dictSources() const
{
    return dictSources_;
}


inline UPtrList<const dictionary>& equationReader::dictSources()
{
    return dictSources_;
}


inline const PtrList<word>& equationReader::dictLookups() const
{
    return dictLookups_;
}


inline const
    UPtrList<const dimensionedScalar>& equationReader::externalDScalars() const
{
    return externalDScalars_;
}


inline UPtrList<const dimensionedScalar>& equationReader::externalDScalars()
{
    return externalDScalars_;
}


inline const UPtrList<const scalar>& equationReader::externalScalars() const
{
    return externalScalars_;
}


inline UPtrList<const scalar>& equationReader::externalScalars()
{
    return externalScalars_;
}


inline const wordList& equationReader::externalScalarNames() const
{
    return externalScalarNames_;
}


inline wordList& equationReader::externalScalarNames()
{
    return externalScalarNames_;
}


inline const
    UPtrList<dimensionedScalar>& equationReader::outputDScalars() const
{
    return outputDScalars_;
}


inline UPtrList<dimensionedScalar>& equationReader::outputDScalars()
{
    return outputDScalars_;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
