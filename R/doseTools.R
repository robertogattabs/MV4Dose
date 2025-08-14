#' class for shaping meshes and other sweets
#'
#' @export
#' @useDynLib MV4Dose
#' @import Morpho Rvcg rgl MV4

doseTools <- function( ) {

  # =============================================================================

  getMeshToPlot <- function( RTDoseMatrix , p.spac ) {

    # RTDoseMatrix e' la matrice 3D che contiene i valori di dose nello spazio, del retto
    # RTDoseMatrix <- r$total.resampled.DoseMSK
    # in p.spac metti il pixel spacing
    # p.spac <- objGL$getPixelSpacing()

    # costruisci gli array con i pixel spacing, per preservare la dimensione nello spazio
    x.arr <- seq(0, (p.spac[1]*(dim(RTDoseMatrix)[1]-1)),by = p.spac[1]  )
    y.arr <- seq(0, (p.spac[2]*(dim(RTDoseMatrix)[2]-1)),by = p.spac[2]  )
    z.arr <- seq(0, (p.spac[3]*(dim(RTDoseMatrix)[3]-1)),by = p.spac[3]  )

    # Metti 0 dove c'erano gli NA
    RTDoseMatrix[which(is.na(RTDoseMatrix),arr.ind = TRUE)] <- -1

    # ora applica una soglia dove il retto comincia a prendere una dose maggiore di 0 Gy
    m.RTDoseMatrix<- contour3d(f = RTDoseMatrix,level = .1, engine="none", x = x.arr, y = y.arr,z = z.arr, plot = FALSE)

    # modifica l'output per renderlo compatibile con il formato mesh
    objS <- services()
    mesh.aaa <- objS$triangle2mesh(x = m.RTDoseMatrix)

    # semplifica la mesh
    mesh.aaa <- vcgClean(mesh = mesh.aaa, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,7,7))
    mesh.aaa <- vcgQEdecim(mesh = mesh.aaa, percent = 0.8)
    mesh.aaa <- vcgSmooth(mesh = mesh.aaa, iteration = 10)

    # ===== COLORAZIONE CON VALORE MASSIMO DEL VICINATO =====

    # Funzione per creare palette dose (rosso scuro -> arancio brillante)
    create_dose_palette <- function() {
      palette_colori <- colorRampPalette(c("#8B0000", "#FF0000", "#FF8C00", "#FFA500"))(100)
      palette_colori <- palette_colori[ length(palette_colori) : 1]
      return( palette_colori )
    }

    # Funzione per trovare il valore massimo del vicinato
    find_max_dose_in_neighborhood <- function(x_coord, y_coord, z_coord,
                                              dose_matrix, x_arr, y_arr, z_arr,
                                              neighborhood_size = 2) {

      # Trova l'indice centrale
      x_idx <- which.min(abs(x_arr - x_coord))
      y_idx <- which.min(abs(y_arr - y_coord))
      z_idx <- which.min(abs(z_arr - z_coord))

      # Calcola le dimensioni della matrice
      dim_x <- dim(dose_matrix)[1]
      dim_y <- dim(dose_matrix)[2]
      dim_z <- dim(dose_matrix)[3]

      # Calcola i range del vicinato
      x_start <- max(1, x_idx - neighborhood_size)
      x_end <- min(dim_x, x_idx + neighborhood_size)
      y_start <- max(1, y_idx - neighborhood_size)
      y_end <- min(dim_y, y_idx + neighborhood_size)
      z_start <- max(1, z_idx - neighborhood_size)
      z_end <- min(dim_z, z_idx + neighborhood_size)

      # Estrai il vicinato e trova il valore massimo
      neighborhood <- dose_matrix[x_start:x_end, y_start:y_end, z_start:z_end]
      max_dose <- max(neighborhood, na.rm = TRUE)

      return(max_dose)
    }

    # Funzione per mappare vertici ai valori di dose usando il vicinato
    map_vertices_to_dose_with_neighborhood <- function(mesh, dose_matrix, x_arr, y_arr, z_arr,
                                                       neighborhood_size = 2) {

      vertices <- mesh$vb[1:3, ]
      n_vertices <- ncol(vertices)
      dose_values <- numeric(n_vertices)

      # Calcola i range delle coordinate
      x_range <- range(x_arr)
      y_range <- range(y_arr)
      z_range <- range(z_arr)

      # Contatori per debug
      inside_range <- 0
      outside_range <- 0
      zero_dose <- 0

      # cat(sprintf("Analizzando %d vertici con vicinato di dimensione %d...\n",
      #             n_vertices, neighborhood_size))

      for (i in 1:n_vertices) {
        x_coord <- vertices[1, i]
        y_coord <- vertices[2, i]
        z_coord <- vertices[3, i]

        # Verifica se il punto è dentro il range della matrice dose
        if (x_coord >= x_range[1] && x_coord <= x_range[2] &&
            y_coord >= y_range[1] && y_coord <= y_range[2] &&
            z_coord >= z_range[1] && z_coord <= z_range[2]) {

          inside_range <- inside_range + 1

          # Trova il valore massimo del vicinato
          dose_value <- find_max_dose_in_neighborhood(x_coord, y_coord, z_coord,
                                                      dose_matrix, x_arr, y_arr, z_arr,
                                                      neighborhood_size)

          dose_values[i] <- dose_value

          if (dose_value == 0) {
            zero_dose <- zero_dose + 1
          }

        } else {
          outside_range <- outside_range + 1
          dose_values[i] <- 0  # Fuori range = dose 0
        }

      }

      return(dose_values)
    }

    # Funzione per colorare la mesh con range ottimizzato
    color_mesh_by_dose_optimized <- function(mesh, dose_values) {

      # Usa i percentili 5% e 95% dei valori > 0 per il range di colorazione
      nonzero_doses <- dose_values[dose_values > 0]

      if (length(nonzero_doses) > 0) {
        min_dose <- quantile(nonzero_doses, 0.05)
        max_dose <- quantile(nonzero_doses, 0.95)
      } else {
        min_dose <- 0
        max_dose <- 1
      }

      # Crea palette
      palette <- create_dose_palette()

      # Normalizza i valori di dose
      if (max_dose > min_dose) {
        normalized_dose <- (dose_values - min_dose) / (max_dose - min_dose)
      } else {
        normalized_dose <- rep(0.5, length(dose_values))
      }

      # Clamp i valori tra 0 e 1
      normalized_dose <- pmax(0, pmin(1, normalized_dose))

      # Mappa ai colori
      color_indices <- round(normalized_dose * (length(palette) - 1)) + 1
      color_indices <- pmin(color_indices, length(palette))

      vertex_colors <- palette[color_indices]
      mesh$material$color <- vertex_colors

      return(list(mesh = mesh,
                  dose_values = dose_values,
                  min_dose = min_dose,
                  max_dose = max_dose,
                  normalized_dose = normalized_dose))
    }

    # ===== APPLICAZIONE COLORAZIONE CON VICINATO =====

    # Testa diverse dimensioni del vicinato
    neighborhood_sizes <- c(1, 2, 3)
    best_result <- NULL
    best_coverage <- 0

    for (size in neighborhood_sizes) {
      # cat(sprintf("\n=== TEST CON VICINATO DI DIMENSIONE %d ===\n", size))

      # Mappa i vertici ai valori di dose usando il vicinato
      vertex_dose_values <- map_vertices_to_dose_with_neighborhood(mesh.aaa, RTDoseMatrix,
                                                                   x.arr, y.arr, z.arr, size)

      # Calcola la copertura (percentuale di vertici con dose > 0)
      coverage <- 100 * sum(vertex_dose_values > 0) / length(vertex_dose_values)

      # cat(sprintf("Copertura dose > 0: %.1f%%\n", coverage))
      # cat(sprintf("Range dose: %.3f - %.3f Gy\n",
      #             min(vertex_dose_values), max(vertex_dose_values)))

      # Scegli il risultato migliore (maggiore copertura)
      if (coverage > best_coverage) {
        best_coverage <- coverage
        best_result <- list(dose_values = vertex_dose_values, size = size)
      }
    }

    # Usa il risultato migliore
    # cat(sprintf("\n=== USANDO VICINATO DI DIMENSIONE %d (copertura %.1f%%) ===\n",
    #             best_result$size, best_coverage))

    # Colora la mesh con range ottimizzato
    colored_mesh_result <- color_mesh_by_dose_optimized(mesh.aaa, best_result$dose_values)

    # Stampa informazioni finali
    # cat(sprintf("Range dose totale: %.3f - %.3f Gy\n",
    #             min(best_result$dose_values), max(best_result$dose_values)))
    # cat(sprintf("Range dose per colorazione: %.3f - %.3f Gy\n",
    #             colored_mesh_result$min_dose, colored_mesh_result$max_dose))

    # Visualizza la mesh colorata
    # shade3d(colored_mesh_result$mesh)

    # Crea legenda
    if (require(fields)) {
      # dev.new()
      dose_range <- seq(colored_mesh_result$min_dose, colored_mesh_result$max_dose, length.out = 100)
      dose_colors <- create_dose_palette()

      legenda <- image.plot(matrix(dose_range, nrow = 1),
                 col = dose_colors,
                 main = sprintf("Legenda Dose (Gy) - Vicinato %dx%dx%d",
                                best_result$size, best_result$size, best_result$size),
                 xlab = "Dose (Gy)",
                 ylab = "",
                 legend.args = list(text = "Dose (Gy)", side = 4, line = 2.5))
    }

    return(  list( "meshToShade" = colored_mesh_result$mesh, "legenda"=legenda ) )

  }

  return( list(
    "getMeshToPlot"=getMeshToPlot
  ))
}


