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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcNone
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    tempSrcField_ = 0.0;
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcStorage
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

#   ifdef FULLDEBUG
    if ((zeroSourceIndex + storageOffset) > maxStoreIndex)
    {
        FatalErrorIn("equationReader::getSouce")
            << "Index " << zeroSourceIndex << " out of bounds (0, "
            << maxStoreIndex - storageOffset << ")"
            << abort(FatalError);
    }
#   endif
    scalarField& returnMe
    (
        storageScalarFields_[zeroSourceIndex + storageOffset]
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcActiveSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    activeSources_[zeroSourceIndex].evaluateScalarField
    (
        tempSrcField_,
        eqOp.componentIndex(),
        geoIndex_
    );
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcEquation
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    dependents_.setSize(dependents_.size() + 1);
    dependents_[dependents_.size() - 1] = equationIndex;

    // Launch the reportEmbeddedDispatchFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedDispatchEnabled;
    //      // or: Info << "Embedded equation dispatch." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedDispatchDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedDispatchFunction_)();

    scalarField result(tempSrcField_.size(), 0.0);
    
    internalEvaluateScalarField
    (
        result,
        zeroSourceIndex,
        maxStoreIndex + 1
    );

    // Launch the reportEmbeddedReturnFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedReturnEnabled;
    //      // or: Info << "Return from equation equation." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedReturnDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedReturnFunction_)();

    tempSrcField_ = result * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcEquationCircRefDetect
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

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
                    operator[](j).name()
                );
                dependencies.append("-->");
            }
            dependencies.append(operator[](i).name());
            FatalErrorIn
            (
                "equationReader::getScalarFieldSrcEquationCircRefDetect"
            )
                << "Circular reference detected when evaluating "
                << "the equation for " << eqn.name()
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

    scalarField result(tempSrcField_.size(), 0.0);
    
    internalEvaluateScalarField
    (
        result,
        zeroSourceIndex,
        maxStoreIndex + 1
    );
    eqOp.assignSourceScalarFunction
    (
        &Foam::equationReader::getScalarSrcEquation
    );
    eqOp.assignSourceScalarFieldFunction
    (
        &Foam::equationReader::getScalarFieldSrcEquation
    );

    if (debug)
    {
        Info << "Returned from embedded equation." << endl;
    }
    tempSrcField_ = result * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcInternalScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    tempSrcField_ = internalScalars_[zeroSourceIndex];
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcDictSourceDScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    dimensionedScalar ds("noSource", dimless, 0.0);
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);
    
    ITstream srcStrm
    (
        dictSources_[zeroSourceIndex].lookup(varName)
    );
    srcStrm >> ds;
    tempSrcField_ = ds.value() * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcDictSourceScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);
    
    scalar returnMe
    (
        readScalar(dictSources_[zeroSourceIndex].lookup(varName))
    );
    tempSrcField_ = returnMe * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcScalarSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    
    tempSrcField_ = scalarSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcScalarFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    
    scalarSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcVectorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    tempSrcField_ = vectorSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcVectorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    vectorSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    tempSrcField_ = tensorSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    
    tensorSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcDiagTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    tempSrcField_ = diagTensorSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcDiagTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    diagTensorSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcSymmTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    tempSrcField_ = symmTensorSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcSymmTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    symmTensorSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcSphericalTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    tempSrcField_ = sphericalTensorSources_.singleValue
    (
        zeroSourceIndex,
        eqOp.componentIndex()
    ) * sign(eqOp.sourceIndex());
    return tempSrcField_;
}


const Foam::scalarField&
    Foam::equationReader::getScalarFieldSrcSphericalTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    sphericalTensorSources_.fullFieldValue
    (
        tempSrcField_,
        zeroSourceIndex,
        eqOp.componentIndex(),
        geoIndex_
    );
    
    tempSrcField_ *= sign(eqOp.sourceIndex());
    return tempSrcField_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
