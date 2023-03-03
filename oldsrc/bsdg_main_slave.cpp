#include "bsdg_mpiflow.h"
#include "mpilead.h"

int
main(int argc, char **argv)
{
  G2_initStandalone(); // TTS 5132362 
  master_t *bsdg = (master_t *)calloc(1, sizeof(master_t));
  bsdg->dmt = (dmt_t *)calloc(1,sizeof(dmt_t));
  bsdg->dmt->v = (struct DMT_MEMORYBUFF_T*)calloc(1, sizeof(struct DMT_MEMORYBUFF_T));
  // bsdg->dbio = (dbio_t *)calloc(1, sizeof(dbio_t));

  // LKL: move this part out of bsdg_main to have separate mpi_init functions for master and slave
  //      To handle mpi_init_thread stuck problem encountered in APAC with timeout
  G2_logInit("bsdg");
  strcpy(bsdg->logName,"bsdg");
  G2_logAddAlias("bsdg",bsdg->logName);

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
  if (provided != MPI_THREAD_SERIALIZED)
      G2_logCritical(bsdg->logName, "MPI doesn't provide MPI_THREAD_SERIALIZED ");

  bsdg_main(argc, argv, bsdg);

  if ( bsdg->l1para.choosemethod == CHOOSEMETHOD_BSDG2D ||
       bsdg->l1para.choosemethod == CHOOSEMETHOD_JOINTSR2D ||
       bsdg->l1para.choosemethod == CHOOSEMETHOD_DAS)
    report_final();
    
  MPI_Finalize();  

  return 0;
}

