#! /usr/bin/env python
import time

from dask.distributed import LocalCluster, Client
print("")

if __name__ == '__main__':

	print("Creating a cluster...")

	worker_opts={'local_dir':'/tmp/dask-worker'}
	#c = LocalCluster(processes=True, n_workers=8, threads_per_worker=1, worker_kwargs=worker_opts)
	#c = LocalCluster(processes=True, n_workers=40, threads_per_worker=1)
	c = LocalCluster(processes=True, n_workers=0, threads_per_worker=1)
	print("Done.")
	print(c)

	print("Scaling...")
	for i in range(40):
		pass
		# Same arguments as the Worker(..) class constructor.
		c.start_worker(ncores=1)
	print("Done.")
	print(c)

	print(c.scheduler)
	#c.scheduler.add_worker(address=c.scheduler.address, local_directory="/tmp/dask-worker")
	#print(c.scheduler)
	#print(c)

	print("Creating a client...")
	client = Client(c)
	print("Done.")

	print(client)


	time.sleep(10)

	print("Closing client...")
	client.close()
	print("Done.")

	print("Closing cluster...")
	c.close()
	print("Done.")
