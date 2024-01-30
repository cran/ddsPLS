# Data Driven Sparse PLS(**ddsPLS**)

**ddsPLS** is a sparse PLS formulation based on soft-thresholding estimations of covariance matrices.

## Installation

There is currently one way to install **ddsPLS**

  * From the under development repository from GitHub thanks to `devtools`

  ```r
  # install.packages("devtools")
  devtools::install_github("hlorenzo/ddsPLS", build_vignettes = TRUE)
  ```
  
Once that package is installed, you can access the vignette using that command.

  ```r
  vignette("ddsPLS")
  ```
  
It is also possible to start a built in applet using 

  ```r
  ddsPLS_App()
  ```

and it should start an interactive interface which should look like

![ddsPLS applet](appCrop.png)

Thanks for using!
