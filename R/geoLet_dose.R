#' class for loading and presenting DICOM Dose data
#'
#' @export
#' @useDynLib MV4Dose
#' @import MV4 XML progress stringr interp 
#' @importFrom mgcv in.out

geoLet_dose <- function(  ) {

  global_objGTL <- c()
  global_lst.Dose <- c()
  global_PSinterpDoseGrid <- c()
  global_cacheDir <- ""
  global_patientFolderName <- ""
  global_resampleDoseVoxelCubeToROI <- list()

  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  # Open a folder and load the content
  addRTDose<-function( pathToOpen ) {

    MM.info <- global_objGTL$getFilesInfo()
    MM.info <- MM.info[which(MM.info[,"kind"] == "RTDoseStorage"),]

    tmpStr <- unlist(strsplit(pathToOpen,"/"))
    global_patientFolderName <<- tmpStr[length(tmpStr)]

    lst.Dose <- list()
    progressivo <- 1
    for( fileName in MM.info[,"fileName"] ) {

      if( file.exists(fileName) == FALSE ) {
        cat("\n Missing ",fileName)
        break;
      }

      righe <- as.numeric(global_objGTL$getTag(tag = "0028,0010" , fileName = fileName ))
      colonne <- as.numeric(global_objGTL$getTag(tag = "0028,0011" , fileName = fileName ))
      PixelSpacing <- global_objGTL$getTag(tag = "0028,0030" , fileName = fileName )
      BitsAllocated <- global_objGTL$getTag(tag = "0028,0100" , fileName = fileName )
      BitsStored <- global_objGTL$getTag(tag = "0028,0101" , fileName = fileName )
      HighBit <- global_objGTL$getTag(tag = "0028,0102" , fileName = fileName )
      PixelRepresentation <- global_objGTL$getTag(tag = "0028,0103" , fileName = fileName )
      DoseUnits <- global_objGTL$getTag(tag = "3004,0002" , fileName = fileName )
      DoseType <- global_objGTL$getTag(tag = "3004,0004" , fileName = fileName )
      GridFrameOffsetVector <- global_objGTL$getTag(tag = "3004,000c" , fileName = fileName )
      FrameOfReferenceUID <- global_objGTL$getTag(tag = "0020,0052" , fileName = fileName )
      SamplesPerPixel <- global_objGTL$getTag(tag = "0028,0002" , fileName = fileName )
      DoseGridScaling <- as.numeric(global_objGTL$getTag(tag = "3004,000e" , fileName = fileName ))
      ImagePositionPatient <- global_objGTL$getTag(tag = "0020,0032" , fileName = fileName )
      ImageOrientationPatient <- global_objGTL$getTag(tag = "0020,0037" , fileName = fileName )
      SOPInstanceUID <- global_objGTL$getTag(tag = "0008,0018" , fileName = fileName )

      numSlices <- length(unlist(strsplit( GridFrameOffsetVector , "\\\\")))
      AbsoluteGridFrameOffsetVector <- as.numeric( unlist(strsplit( GridFrameOffsetVector , "\\\\")) )
      AbsoluteGridFrameOffsetVector <- AbsoluteGridFrameOffsetVector + as.numeric(unlist(strsplit( ImagePositionPatient , "\\\\"))[3])

      if( SamplesPerPixel != "1") { stop(" SamplesPerPixel is expected to be '1' in the current version of moddicom") }
      if( DoseType != "PHYSICAL") { stop(" SamplesPerPixel is expected to be 'PHYSICAL' in the current version of moddicom") }
      if( DoseUnits != "GY") { stop(" SamplesPerPixel is expected to be 'GY' in the current version of moddicom") }
      if( PixelRepresentation != "0") { stop(" PixelRepresentation is expected to be '0' in the current version of moddicom") }
      if( HighBit != "31") { stop(" HighBit is expected to be '31' in the current version of moddicom") }
      if( BitsStored != "32") { stop(" BitsStored is expected to be '32' in the current version of moddicom") }
      if( BitsAllocated != "32") { stop(" HighBit is expected to be '32' in the current version of moddicom") }

      # load the RTDose
      cat("\n Loading: ",fileName)
      rn <- readBin(con = fileName, what = "integer", size = 1, endian = "little",
                    n = file.size(fileName), signed = FALSE)
      rn <- rn[ length(rn):1 ]

      # assign the byte to the proper Dose bin
      cat("\n\tProcessing: ",fileName)
      aaa <- unlist(lapply( seq(1, ( righe * colonne * numSlices * 4) , by = 4) , function( pos ) {
        rn[ (pos+3) ] + rn[ (pos+2) ] * 2^8 + rn[ (pos+1) ] * 2^16 + rn[ (pos+0) ] * 2^24
      }))
      rn <- aaa

      # prepare the matrix
      oppa <- rn[ 1:( righe * colonne * numSlices)   ]
      oppa <- oppa[ length(oppa):1 ]
      matRN <- array(0,c(righe,colonne,numSlices))

      # fill the matrix
      # ct<-1
      # for( z in seq(1,numSlices)) {
      #   for(x in seq(1,righe)) {
      #     for(y in seq(1,colonne)) {
      #       matRN[x,colonne-y,z]<-oppa[ct]
      #       ct<-ct+1
      #     }
      #   }
      # }
      ct<-1
      t <- lapply( seq(1,numSlices), function(z){
        t <- lapply( seq(1,righe), function(x){
          t <- lapply( seq(1,colonne), function(y){
            matRN[x,colonne-y,z] <<- oppa[ct]
            ct <<- ct+1
            } )
          } )
      } )

      # Ribalta il voxel cube cosi' da averlo coerente con il VC della TAC
      new.matRN <- array(0,dim =  dim(matRN)[c(2,1,3)] )
      tmp <- lapply( 1:(dim(matRN)[3]), function(z){  new.matRN[,,z] <<- t(matRN[dim(matRN)[1]:1,,z]) })
      matRN <- new.matRN

      # Non dimentricare il gridScaling
      matRN <- matRN * DoseGridScaling

      # chiave <- as.character(progressito)
      chiave <- SOPInstanceUID

      lst.Dose[[ chiave ]] <- list()
      lst.Dose[[ chiave  ]]$FrameOfReferenceUID <- FrameOfReferenceUID
      lst.Dose[[ chiave  ]]$Dose <- matRN
      lst.Dose[[ chiave  ]]$GridFrameOffsetVector <- GridFrameOffsetVector
      lst.Dose[[ chiave  ]]$AbsoluteGridFrameOffsetVector <- AbsoluteGridFrameOffsetVector
      lst.Dose[[ chiave  ]]$ImagePositionPatient <- unlist(strsplit(ImagePositionPatient, "\\\\"))
      lst.Dose[[ chiave  ]]$ImageOrientationPatient <- unlist(strsplit(ImageOrientationPatient, "\\\\"))
      lst.Dose[[ chiave  ]]$PixelSpacing <- unlist(strsplit(PixelSpacing, "\\\\"))

      progressivo <- progressivo + 1
    }

    global_lst.Dose <<- lst.Dose
  }


getROIDose <- function( ROIName , CT.SeriesInstanceUID = NA , cropIt = TRUE ) {

    arr.RD.SOPInstanceUID <- get.RD.SOPInstanceUID();
    dataStorage  <- global_objGTL$getDataStorage();

    # per ogni SOPInstanceUID di dose, calcola la dose prevista
    lst.stacked.MMask <- list()
    total.resampled.DoseVC <- c()
    for( SOPInstanceUID in arr.RD.SOPInstanceUID  ) {
      doseStorage <- global_lst.Dose[[SOPInstanceUID]]

      # if( global_cacheDir == "" ) stop("\n cacheDir not set ")

      # global_cacheDir.fullName <- paste(c(global_cacheDir,"/",global_patientFolderName),collapse='')
      # #
      # if( dir.exists(global_cacheDir.fullName) == FALSE ) {
      #   dir.create( global_cacheDir.fullName )
      #   if( dir.exists(global_cacheDir.fullName) == FALSE ) {
      #     stop("Dir not created: Error ##h98hh")
      #   }
      # }

      # Prepara i dati per avere i valori di dose interpolati sui punti della TC

      CT.VC <- global_objGTL$getImageVoxelCube(SeriesInstanceUID = CT.SeriesInstanceUID)
      CT.filesInfo <- global_objGTL$getFilesInfo()
      CT.filesInfo <- CT.filesInfo[ which(CT.filesInfo[,"SeriesInstanceUID"] == CT.SeriesInstanceUID) , ]
      CT.filesInfo <- CT.filesInfo[order(as.numeric(CT.filesInfo[,"IPP.z"])),]
      IPP <- as.numeric(doseStorage$ImagePositionPatient)
      IOP <- as.numeric(doseStorage$ImageOrientationPatient)
      x.start <- IPP[1];  y.start <- IPP[2]

      arr.x.CT <- seq(0:(dim(CT.VC)[1]-1)) * (as.numeric(CT.filesInfo[1,"p.x"]))
      arr.y.CT <- seq(0:(dim(CT.VC)[2]-1)) * (as.numeric(CT.filesInfo[1,"p.y"]))
      arr.x.CT <- arr.x.CT + as.numeric(CT.filesInfo[1,"IPP.x"])
      arr.y.CT <- arr.y.CT + as.numeric(CT.filesInfo[1,"IPP.y"])
      arr.z.CT <- sort(as.numeric(CT.filesInfo[,"IPP.z"]))

      arr.x.DS <- seq(0:(dim(doseStorage$Dose)[1]-1)) * as.numeric(doseStorage$PixelSpacing)[1]
      arr.y.DS <- seq(0:(dim(doseStorage$Dose)[2]-1)) * as.numeric(doseStorage$PixelSpacing)[2]
      arr.x.DS <- arr.x.DS + x.start
      arr.y.DS <- arr.y.DS + y.start
      arr.z.DS <- sort(doseStorage$AbsoluteGridFrameOffsetVector)

      # Invoca l'interpolazione trilineare in ANSI C
      a <- interpola.C(arr.x.CT = arr.x.CT, arr.y.CT = arr.y.CT, arr.z.CT = arr.z.CT,
                       arr.x.DS = arr.x.DS, arr.y.DS = arr.y.DS, arr.z.DS = arr.z.DS,
                       Matrice3D = doseStorage$Dose)

      quali <- which(is.na(a),arr.ind = T)
      if( length(quali) > 0 ) {
        a[quali] <-  0
      }

      # incrementa la matrice 3D con i totali di dose ricampionati
      if( SOPInstanceUID == arr.RD.SOPInstanceUID[1] ) {
        total.resampled.DoseVC  <- a
      } else {
        total.resampled.DoseVC <- total.resampled.DoseVC + a
      }

    }

    # ok, ora mi aspetto che le slice siano complanari alle ROI: inizio ad affrontare le polyline
    # MMM <- total.resampled.DoseVC
    dim.TRDVC <- dim(total.resampled.DoseVC)
    arr.z.ROI <- arr.z.CT
    stacked.MMask <- array( 0 ,dim = dim.TRDVC )
    for(SOPQualcosa in names(dataStorage$structures[[ROIName]]) ) {

      quante.polyline <- length(dataStorage$structures[[ROIName]][[SOPQualcosa]])
      partial.stacked.MMask <- array( 0 ,dim = c(dim.TRDVC[1],dim.TRDVC[2],quante.polyline)  )
      contatore <- 1
      for( n.polyline in 1:quante.polyline ) {

        # Per ora prendi solo la prima polyline, poi questo '1' andra' in un loop'
        ROI.point.table <- dataStorage$structures[[ROIName]][[SOPQualcosa]][[n.polyline]]

        # cerca la slice corrispondente
        current.z  <- ROI.point.table[1,3]
        sliceToConsider <- which(arr.z.ROI == current.z)

        # estrai i dati per il point In Polygon
        arr.x.coor <- (0:(dim.TRDVC[1]-1)) * as.numeric(doseStorage$PixelSpacing[1]) + as.numeric( doseStorage$ImagePositionPatient[1] )
        arr.y.coor <- (0:(dim.TRDVC[2]-1)) * as.numeric(doseStorage$PixelSpacing[2]) + as.numeric( doseStorage$ImagePositionPatient[2] )

        PSP <- unlist(strsplit(global_objGTL$getTag("0028,0030",SeriesInstanceUID = CT.SeriesInstanceUID),"\\\\"))
        IPP <- unlist(strsplit(global_objGTL$getTag("0020,0032",SeriesInstanceUID = CT.SeriesInstanceUID),"\\\\"))
        arr.x.coor <- (0:(dim.TRDVC[1]-1)) * as.numeric(PSP[1]) + as.numeric( IPP[1] )
        arr.y.coor <- (0:(dim.TRDVC[2]-1)) * as.numeric(PSP[2]) + as.numeric( IPP[2] )

        points2Test <- as.matrix(expand.grid( arr.x.coor , arr.y.coor ))
        ROIPoints <- ROI.point.table[,1:2]

        inside <- which(in.out(ROIPoints,as.matrix(points2Test)))

        # Costruisci la matrice maschera
        MMask <- matrix( 0, nrow = dim.TRDVC[1], ncol = dim.TRDVC[2] )
        rownames(MMask) <- as.character(arr.x.coor)
        colnames(MMask) <- as.character(arr.y.coor)

        # Se ci sono punti interni, aggiungili
        if( length(inside) > 0 ) {
          tmp <- apply(X = points2Test[inside,],MARGIN = 1,FUN = function(i){  MMask[ as.character(i[1]) , as.character(i[2])] <<- 1  })
        }

        MMask  <- (MMask[dim(MMask)[1]:1,dim(MMask)[2]:1])
        # stacked.MMask[,,sliceToConsider] <- MMask
        partial.stacked.MMask[,,n.polyline] <- MMask

        contatore <- contatore + 1
      }
      # componi le partial.stacked.MMask per costruire la stacked.MMask
      # if( n.polyline > 1) {
      #   a <- apply(partial.stacked.MMask,MARGIN = c(1,2) ,sum)
      #   tmp  <- which( (a %% 2) ==0 ,arr.ind = T)
      #   if( length(tmp) > 0 ) { a[tmp] <- NA }
      #   tmp <- which(!is.na(a),arr.ind = T)
      #   if( length(tmp) > 0 ) { a[tmp] <- 1  }
      #
      #   mtr.new.contribute <- apply(partial.stacked.MMask, MARGIN = c(1,2),sum)
      #
      #   # Se voglio tenerli tutti:
      #   mtr.new.contribute[which( mtr.new.contribute > 0, arr.ind = T)] <- 1
      #   # Se voglio tenere solo i dispari
      #   # mtr.new.contribute[which( (mtr.new.contribute %% 2) == 1 , arr.ind = T)] <- 1
      #   # mtr.new.contribute[which( (mtr.new.contribute %% 2) != 1 , arr.ind = T)] <- 0
      #
      #   stacked.MMask[,,sliceToConsider] <- stacked.MMask[,,sliceToConsider] + mtr.new.contribute
      # } else{
      #   stacked.MMask[,,sliceToConsider] <- partial.stacked.MMask[,,1]
      # }
# browser()
      mtr.new.contribute <- apply(partial.stacked.MMask, MARGIN = c(1,2),sum)
      # Prova a tenere solo i dispari
      mtr.new.contribute[which( (mtr.new.contribute %% 2) == 1 , arr.ind = T)] <- 1
      mtr.new.contribute[which( (mtr.new.contribute %% 2) != 1 , arr.ind = T)] <- 0

      # se non riesci (che si cancella tutto), prova a tenere tutti
      # TODO: Testare se ha senso questa logica: e' necessario avere piu' RTDose e RTStruct con delle roi cave
      if( sum(mtr.new.contribute) == 0 ) {
        mtr.new.contribute <- apply(partial.stacked.MMask, MARGIN = c(1,2),sum)
        mtr.new.contribute[which( mtr.new.contribute > 0, arr.ind = T)] <- 1
      }

      # aggiungi alle slice fatte
      stacked.MMask[,,sliceToConsider] <- mtr.new.contribute

    }

    # browser()

    total.resampled.DoseMSK <- total.resampled.DoseVC * stacked.MMask
    total.resampled.DoseMSK[which( total.resampled.DoseMSK == 0, arr.ind = TRUE)] <- NA
    total.resampled.DoseMSK <- total.resampled.DoseMSK[ dim(total.resampled.DoseMSK)[1]:1,, ]

    BIDIDOSONE <- total.resampled.DoseMSK[!is.na(total.resampled.DoseMSK)]
    dose.max <- max(BIDIDOSONE)
    x.hist <- seq(0,dose.max,by=(dose.max/100))
    y.hist <- unlist(lapply( x.hist, function(soglia){  length(which(BIDIDOSONE>=soglia))/length(BIDIDOSONE)    }   ))

    original.CT <- global_resampleDoseVoxelCubeToROI$CT.MMM
    arr.z.coords <- global_resampleDoseVoxelCubeToROI$arr.z.ROI


    return( list("lst.stacked.MMask" = lst.stacked.MMask,
                 "x.hist"=x.hist ,"y.hist"=y.hist,
                 "total.resampled.DoseMSK"=total.resampled.DoseMSK,
                 "original.CT" = original.CT,
                 "arr.z.coords" = arr.z.coords))
  }



  resampleDoseVoxelCubeToROI <- function( arr.RD.SOPInstanceUID , ROIName, CT.SeriesInstanceUID = NA ) {
    # prendi delle variabili di uso generale
    dataStorage <- global_objGTL$getDataStorage()

    # Prendi le coordinate z della ROI di interesse
    arr.z <- unlist(lapply( names(dataStorage$structures[[ROIName]]), function(chiave){
      quali.z <- lapply( 1:length(dataStorage$structures[[ROIName]][[chiave]]) , function(i){
        z.pos <- unique(dataStorage$structures[[ROIName]][[chiave]][[i]][,3])
        if( length(z.pos) > 1 ) { browser(); stop("ERROR: only axial ROIs are allowed")}
        return( z.pos )
      } )
      quali.z <- unique(unlist(quali.z))
      if( length(quali.z) > 1 ){ browser(); stop("ERROR: only axial ROIs are allowed") }
      return( quali.z );
    } ))

    # browser()

    # metti i dati in una struttura tabellare: l'ssunto di base e' che le ROI siano state costruite sulle CT
    mtr.roi.info <- cbind( names(dataStorage$structures[[ROIName]]), arr.z)
    if( !is.na(CT.SeriesInstanceUID) ) {
      CT.filesInfo <- global_objGTL$getFilesInfo()
      CT.filesInfo <- CT.filesInfo[ which(CT.filesInfo[,"SeriesInstanceUID"] == CT.SeriesInstanceUID) , ]
      CT.filesInfo <- CT.filesInfo[order(as.numeric(CT.filesInfo[,"IPP.z"])),]
      CT.VC <- global_objGTL$getImageVoxelCube(SeriesInstanceUID = CT.SeriesInstanceUID)
      # CT.MMM <- array( 0,dim = c( dim(CT.VC[,,1]), nrow(mtr.roi.info)  ))
      CT.MMM <- array( 0, dim = dim(CT.VC) )
      if( dim(CT.VC)[3] != nrow(CT.filesInfo)  ) {
        CT.error <- true
        cat("ERROR:   dim(CT.VC)[3] != nrow(CT.filesInfo)")
        stop();
      }
    }



    # prendi la dose di quella specifica SOPInstanceUID
    SOPInstanceUID <- arr.RD.SOPInstanceUID[[1]]
    doseStorage <- global_lst.Dose[[SOPInstanceUID]]

    # ora per ogni polylinec cerca quali slide di DOSE devono essere considerate, superiormente e inferiormente
    arr.z.ROI <- c()
    # MMM <- array( 0,dim = c( dim(doseStorage$Dose[,,1]), nrow(mtr.roi.info)  ))
    # MMM.interp <- array( 0,dim = c( dim(CT.VC[,,1]), nrow(mtr.roi.info)  ))
    MMM <- array( 0,dim = dim(doseStorage$Dose) )
    MMM.interp <- array( 0, dim = dim(CT.VC) )

    # for( riga in 1:nrow(mtr.roi.info)) {
    num.slice.CT <- dim(CT.VC)[3]
    for( riga in 1:num.slice.CT ) {
      # -------------------------------------------------------------
      # PER LA MATRICE DI DOSE
      # -------------------------------------------------------------
      z <- as.numeric(mtr.roi.info[riga,2]  )
      x <- doseStorage$AbsoluteGridFrameOffsetVector - z

      # Trova gli indici dove c'è un cambio di segno
      pos_change_indices <- which( (x>0) != (x[1]>0) )[1]-1

      # SE e' fuori campo, cattura l'eccezione
      if( is.na(pos_change_indices) ) {
        browser()
        stop()
      }

      IPP <- as.numeric(doseStorage$ImagePositionPatient)
      IOP <- as.numeric(doseStorage$ImageOrientationPatient)
      x.start <- IPP[1];  y.start <- IPP[2]
      browser()
      # In un caso non devo interpolare
      if( x[pos_change_indices] == 0 ) {
        MMM[,,riga] <- doseStorage$Dose[,,pos_change_indices]
        arr.z.ROI <- c( arr.z.ROI   ,  z)
      } else {
        # In un altro, be', si'
        z.immagine.over <- doseStorage$AbsoluteGridFrameOffsetVector[pos_change_indices]
        z.immagine.under <-doseStorage$AbsoluteGridFrameOffsetVector[pos_change_indices+1]
        z.ROI <- z

        slice.immagine.over <- which( doseStorage$AbsoluteGridFrameOffsetVector  == z.immagine.over)
        slice.immagine.under<- which( doseStorage$AbsoluteGridFrameOffsetVector  == z.immagine.under)

        # IPP <- as.numeric(doseStorage$ImagePositionPatient)
        # IOP <- as.numeric(doseStorage$ImageOrientationPatient)
        #
        # x.start <- IPP[1];  y.start <- IPP[2]
        # browser()
        arr.x.img <- seq(0:(dim(doseStorage$Dose)[1]-1)) * as.numeric(doseStorage$PixelSpacing)[1]
        arr.y.img <- seq(0:(dim(doseStorage$Dose)[2]-1)) * as.numeric(doseStorage$PixelSpacing)[2]
        arr.x.img <- arr.x.img + x.start
        arr.y.img <- arr.y.img + y.start

        # Visto che i punti sono esattamente sottostanti, non serve che scomodi una trilineare

        DX <- (z.immagine.over - z.immagine.under)
        mtr.m <- (doseStorage$Dose[,,slice.immagine.over] - doseStorage$Dose[,,slice.immagine.under])/DX

        x <- z.ROI - z.immagine.under
        mtr.int <- (mtr.m * x )+ doseStorage$Dose[,,slice.immagine.under]

        # image( doseStorage$Dose[,,slice.immagine.over]  )
        # image( mtr.int  )
        # image( doseStorage$Dose[,,slice.immagine.under]  )

        arr.z.ROI <- c( arr.z.ROI   ,  z)
        MMM[,,riga] <- mtr.int
      }
      # -------------------------------------------------------------
      # PER LA CT
      # -------------------------------------------------------------
      # browser()
      if( !is.na(CT.SeriesInstanceUID) ) {

        CT.x <- as.numeric(CT.filesInfo[,"IPP.z"])  - z
        CT.pos_change_indices <- which( (CT.x>0) != (CT.x[1]>0) )[1]-1

        # In un caso non devo interpolare
        if( CT.x[CT.pos_change_indices] == 0 ) {
          # CT.MMM[,,riga] <- CT.VC[,,pos_change_indices]
          CT.MMM[,,riga] <- CT.VC[,,CT.pos_change_indices]
        } else {
          browser()
        }

        # Qui calcola la matrice di dose interpolata sulla CT

        arr.x.CT <- seq(0:(dim(CT.VC)[1]-1)) * (as.numeric(CT.filesInfo[1,"p.x"]))
        arr.y.CT <- seq(0:(dim(CT.VC)[2]-1)) * (as.numeric(CT.filesInfo[1,"p.y"]))
        arr.x.CT <- arr.x.CT + as.numeric(CT.filesInfo[1,"IPP.x"])
        arr.y.CT <- arr.y.CT + as.numeric(CT.filesInfo[1,"IPP.y"])
        # browser()
        arr.x.img <- seq(0:(dim(doseStorage$Dose)[1]-1)) * as.numeric(doseStorage$PixelSpacing)[1]
        arr.y.img <- seq(0:(dim(doseStorage$Dose)[2]-1)) * as.numeric(doseStorage$PixelSpacing)[2]
        arr.x.img <- arr.x.img + x.start
        arr.y.img <- arr.y.img + y.start
        x.arr.Dose <- arr.x.img
        y.arr.Dose <- arr.y.img

        tmp.mtr <- bilinear.grid(x = x.arr.Dose, y = y.arr.Dose, z = MMM[,,riga],
                                 xlim = range(arr.x.CT) , ylim = range(arr.y.CT),
                                 nx = length(arr.x.CT),ny = length(arr.y.CT)  )$z

        MMM.interp[,,riga] <- tmp.mtr
      }
      # pb$tick()
    }
    # per il momento ritornane solo una, piu' avanti potresti restituirne una per ogni matricedi gruppo di polyline'
    # se decidero' di splittarle'

    global_resampleDoseVoxelCubeToROI <<-
      list( "MMMs" = list("1" = MMM), "arr.z.ROI"=list("1"=arr.z.ROI),
            "MMMs.interp" = list("1" = MMM.interp),
            "CT.MMM" = CT.MMM
      )

  }


  set_objGTL <- function( objGTL ) {
    global_objGTL <<- objGTL
  }
  set_resampleDoseVoxelCubeToROI <- function( given_resampleDoseVoxelCubeToROI ) {
    global_resampleDoseVoxelCubeToROI <<- given_resampleDoseVoxelCubeToROI
  }
  getDateStorage <- function() {
    return( list(
      "global_lst.Dose"=global_lst.Dose
    ))
  }
  get.RD.SOPInstanceUID <- function() {
    return( names(global_lst.Dose))
  }
  #=================================================================================
  # Constructor
  #=================================================================================




  interpola.C <- function( arr.x.CT , arr.y.CT , arr.z.CT,
                           arr.x.DS , arr.y.DS, arr.z.DS,
                           Matrice3D ) {

    # interpola.C <- function( xPos_CT , yPos_CT , zPos_CT,
    #                          xPos_DS , yPos_DS, zPos_DS,
    #                          Matrice3D ) {


    # Matrice3D <- doseStorage$Dose
    newXps <- as.vector(aperm(Matrice3D, c(3, 2, 1)))

    n_out <- length(arr.x.CT) * length(arr.y.CT) * length(arr.z.CT)

    result <- .C("interpola_suAltraMatricePunti",
                 as.double(arr.x.CT), as.double(arr.y.CT), as.double(arr.z.CT),
                 as.double(arr.x.DS), as.double(arr.y.DS), as.double(arr.z.DS),
                 as.integer(length(arr.x.CT)), as.integer(length(arr.y.CT)), as.integer(length(arr.z.CT)),
                 as.integer(length(arr.x.DS)), as.integer(length(arr.y.DS)), as.integer(length(arr.z.DS)),
                 as.double(newXps),
                 double(n_out)
    )

    a <- result[[14]]

    img <- aperm(array(a, dim = c(length(arr.z.CT), length(arr.y.CT), length(arr.x.CT))), c(3, 2, 1))

    return(img)




    # newXps <- as.vector(Matrice3D)
    #
    # n_out <- length(xPos_CT) * length(yPos_CT) * length(zPos_CT)
    #
    # result <- .C("interpola_suAltraMatricePunti",
    #              as.double(xPos_CT), as.double(yPos_CT), as.double(zPos_CT),
    #              as.double(xPos_DS), as.double(yPos_DS), as.double(zPos_DS),
    #              as.integer(length(xPos_CT)), as.integer(length(yPos_CT)), as.integer(length(zPos_CT)),
    #              as.integer(length(xPos_DS)), as.integer(length(yPos_DS)), as.integer(length(zPos_DS)),
    #              as.double(newXps),
    #              double(n_out)
    # )[[14]]
    #
    # browser()
    #
    # a <- result
    #
    # objService <- services()
    # matrice <- array( result , dim=c( length(yPos_CT) , length(xPos_CT), length(zPos_CT) ))
    # for ( i in seq(1,dim(matrice)[3] )) {
    #   matrice[,,i]<-t(objService$rotateMatrix(matrice[,,i],rotations=3))
    # }

    # browser()
    #
    # a <- 2
    # a <- 2
    # a <- 2
    # a <- 2



  }

  # sommalo.C <- function() {
  #   a <- 1.2
  #   b <- 6.5
  #   ccc <- 0
  #   .C("sommalo",  as.double(a),as.double(b), as.double(ccc))
  #   browser()
  #   a <- 4
  #   a <- 7
  # }


  # test_organ_detection <- function() {
  # 
  #   CT.VC <- global_objGTL$getImageVoxelCube()
  # 
  #   lst.res <- list()
  #   for( slice in 1:dim(CT.VC)[3] ) {
  #     # Array per risultati: [polmoni, cuore, encefalo, vescica]
  #     results <- integer(4)
  # 
  #     # Chiama la funzione C
  #     result <- .C("organ_detection_wrapper",
  #                  as.double(slice),
  #                  as.integer(dim(CT.VC)[1]),
  #                  as.integer(dim(CT.VC)[2]),
  #                  results = as.integer(results))
  # 
  #     # Estrai i risultati
  #     results <- result$results
  # 
  #     # Nomi degli organi
  #     # organ_names <- c("Polmoni", "Cuore", "Encefalo", "Vescica")
  # 
  #     # Stampa risultati
  #     # cat("=== RISULTATI RICONOSCIMENTO ORGANI ===\n")
  #     # for (i in 1:4) {
  #     #   status <- if (results[i]) "PRESENTE" else "ASSENTE"
  #     #   cat(sprintf("%s: %s\n", organ_names[i], status))
  #     # }
  #     lst.res[[ as.character(slice) ]] <- results
  # 
  #   }
  # 
  # 
  # 
  #   return(lst.res)
  # }

  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function(  ) {
    global_objGTL <<- c()
    global_lst.Dose <<- c()
    global_PSinterpDoseGrid <<- c( 0.75, 0.75, 1 )
    global_resampleDoseVoxelCubeToROI <<- list()
    global_cacheDir <<- ""
  }
  constructor(  )
  return( list(
    "addRTDose"=addRTDose,
    "getROIDose"=getROIDose,
    "set_objGTL"=set_objGTL,
    "getDateStorage"=getDateStorage,
    "get.RD.SOPInstanceUID"=get.RD.SOPInstanceUID,
    "resampleDoseVoxelCubeToROI"=resampleDoseVoxelCubeToROI,
    "set_resampleDoseVoxelCubeToROI"=set_resampleDoseVoxelCubeToROI
    # "test_organ_detection"=test_organ_detection
  ))


}

# $env:RTOOLS="C:\Rtools";
# $env:PATH="$env:RTOOLS\bin;$env:RTOOLS\mingw64\bin;$env:PATH"
# cmd /c "R CMD SHLIB base.c"

#
#
# library(MV4)
# library(mgcv)
# library(progress)
# library(interp)
# library(scales)
# library(stringr)
#
# lst.meshes.fileName <- "./tmp/lst.lst.meshes.RData"
# arr.patientFolderName  <- c("AA08","AA16","FV07","FA30","JR14","JR23","JS25","JU02")
# arr.patientFolderName  <- c("JR23","AA08","AA16","FV07","FA30","JR14","JR23","MA06","MA12","MA29","MD01")
#
# if( file.exists(lst.meshes.fileName)) {
#   load( lst.meshes.fileName )
# } else {
#   lst.meshes <- list()
# }
#
# objS <- services()
# for( patientFolderName in arr.patientFolderName )  {
#
#   if( !(patientFolderName %in% names(lst.meshes)) ) {
#     patientGeoLetFileName <- paste(c("./tmp/",patientFolderName,".RData"),collapse = '')
#
#     if( file.exists(patientGeoLetFileName) ) {
#       load( patientGeoLetFileName  )
#     } else {
#       patientDICOMFileName <- paste(c("C://projects/images/dosomicsNavarra/",patientFolderName),collapse = '')
#       objGL <- geoLet()
#       objGL$openDICOMFolder(pathToOpen = patientDICOMFileName )
#       save( objGL , file = patientGeoLetFileName )
#
#       tmp <- file.remove( list.files(patientDICOMFileName,pattern = ".xml",full.names = T) )
#       tmp <- file.remove( list.files(patientDICOMFileName,pattern = ".raw",full.names = T) )
#     }
#
#
#     GLTD <- geoLet_dose( cacheDir = "./tmp/" )
#     GLTD$set_objGTL( objGL )
#     GLTD$addRTDose(pathToOpen = paste(c("c:\\projects/images/dosomicsNavarra/",patientFolderName,"/"),collapse = '')  )
#
#
#     CT.SeriesInstanceUID <- objGL$get.CT.SeriesInstanceUID()[1]
#
#     ROIName <- "GTVn1"
#     # get the 3D dose distribution
#     r <- GLTD$getROIDose( ROIName = ROIName,CT.SeriesInstanceUID =  CT.SeriesInstanceUID, cropIt = FALSE)
#
#     stop()
#
#     # get the mesh
#     aaa <- r$cropped.DOSE
#     p.spac <- objGL$getPixelSpacing()
#     x.arr <- seq(0, (p.spac[1]*(dim(aaa)[1]-1)),by = p.spac[1]  )
#     y.arr <- seq(0, (p.spac[2]*(dim(aaa)[2]-1)),by = p.spac[2]  )
#     z.arr <- seq(0, (p.spac[3]*(dim(aaa)[3]-1)),by = p.spac[3]  )
#     aaa[which(is.na(aaa),arr.ind = TRUE)] <- 0
#
#     m.aaa <- contour3d(f = aaa,level = .1, engine="none", x = x.arr, y = y.arr,z = z.arr)
#     mesh.aaa <- objS$triangle2mesh(x = m.aaa)
#     mesh.aaa<-vcgClean(mesh = mesh.aaa, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,7,7))
#     mesh.aaa <- vcgQEdecim(mesh = mesh.aaa, percent = 0.8)
#     mesh.aaa <- vcgSmooth(mesh = mesh.aaa, iteration = 10)
#
#
#
#     lst.meshes[[ patientFolderName ]]$ROI <- list()
#     lst.meshes[[ patientFolderName ]]$ROIS$Rectum.volume <- r
#     lst.meshes[[ patientFolderName ]]$ROIS$Rectum.mesh <- mesh.aaa
#
#     save( lst.meshes , file = lst.meshes.fileName )
#   }
#
# }
#
#
