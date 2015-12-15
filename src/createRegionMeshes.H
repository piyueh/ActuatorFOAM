const word          fluidName(rp["fluid"][0]);
const word          solidName(rp["solid"][0]);


Info << "Create default mesh for time = "                                             
     << runTime.timeName() << nl << endl;                         
                                                                                  
fvMesh defaultMesh                                                            
(                                                                            
    IOobject                                                           
    (                                                                        
        fvMesh::defaultRegion,                                         
        runTime.timeName(),                                                  
        runTime,                                                             
        IOobject::MUST_READ                                            
    )                                                                        
);


Info << "Create fluid mesh for region " 
     << fluidName << " for time = "
     << runTime.timeName() << nl << endl;

fvMesh fluidMesh
(
    IOobject
    (
        fluidName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);


Info << "Create solid mesh for region " 
     << solidName << " for time = "
     << runTime.timeName() << nl << endl;

fvMesh solidMesh
(
    IOobject
    (
        solidName,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);

