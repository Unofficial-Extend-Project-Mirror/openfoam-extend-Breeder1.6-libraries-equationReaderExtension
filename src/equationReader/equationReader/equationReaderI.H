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

inline const Foam::PtrList<scalar>&
    Foam::equationReader::internalScalars() const
{
    return internalScalars_;
}


inline const Foam::UPtrList<const dictionary>&
    Foam::equationReader::dictSources() const
{
    return dictSources_;
}


inline const Foam::PtrList<word>& Foam::equationReader::dictLookups() const
{
    return dictLookups_;
}


inline const Foam::UPtrList<const equationVariable>&
    Foam::equationReader::activeSources() const
{
    return activeSources_;
}


inline const Foam::wordList& Foam::equationReader::activeSourceNames() const
{
    return activeSourceNames_;
}


inline const Foam::PtrList<dimensionSet>&
    Foam::equationReader::activeSourceDimensions() const
{
    return activeSourceDimensions_;
}


inline const Foam::equationSource<scalar>&
    Foam::equationReader::scalarSources() const
{
    return scalarSources_;
}


inline Foam::equationSource<scalar>& Foam::equationReader::scalarSources()
{
    return scalarSources_;
}


inline const Foam::equationSource<vector>&
    Foam::equationReader::vectorSources() const
{
    return vectorSources_;
}


inline Foam::equationSource<vector>& Foam::equationReader::vectorSources()
{
    return vectorSources_;
}


inline const Foam::equationSource<tensor>&
    Foam::equationReader::tensorSources() const
{
    return tensorSources_;
}


inline Foam::equationSource<tensor>& Foam::equationReader::tensorSources()
{
    return tensorSources_;
}


inline const Foam::equationSource<diagTensor>&
    Foam::equationReader::diagTensorSources() const
{
    return diagTensorSources_;
}


inline Foam::equationSource<diagTensor>&
    Foam::equationReader::diagTensorSources()
{
    return diagTensorSources_;
}


inline const Foam::equationSource<symmTensor>&
    Foam::equationReader::symmTensorSources() const
{
    return symmTensorSources_;
}


inline Foam::equationSource<symmTensor>&
    Foam::equationReader::symmTensorSources()
{
    return symmTensorSources_;
}


inline const Foam::equationSource<sphericalTensor>&
    Foam::equationReader::sphericalTensorSources() const
{
    return sphericalTensorSources_;
}


inline Foam::equationSource<sphericalTensor>&
    Foam::equationReader::sphericalTensorSources()
{
    return sphericalTensorSources_;
}


/*inline const Foam::label& Foam::equationReader::geoIndex() const
{
    return geoIndex_;
}


inline void Foam::equationReader::setGeoIndex(label newIndex)
{
    geoIndex_ = newIndex;
}


inline const Foam::label& Foam::equationReader::cellIndex() const
{
    return cellIndex_;
}


inline void Foam::equationReader::setCellIndex(label newIndex)
{
    cellIndex_ = newIndex;
}*/

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
