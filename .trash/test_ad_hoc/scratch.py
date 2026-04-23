import asyncio
import time

async def main():
    q = asyncio.Queue(maxsize=3)
    def prod(loop, queue):
        print("prod start")
        for i in range(5):
            print("putting", i)
            asyncio.run_coroutine_threadsafe(queue.put(i), loop).result()
            print("put", i)
        asyncio.run_coroutine_threadsafe(queue.put(None), loop).result()
        print("prod end")

    loop = asyncio.get_running_loop()
    producer_task = asyncio.to_thread(prod, loop, q)

    while True:
        print("waiting for item")
        item = await q.get()
        print("got", item)
        if item is None:
            break
        await asyncio.sleep(0.1)

    print("waiting for producer task")
    await producer_task
    print("done")

asyncio.run(main())
