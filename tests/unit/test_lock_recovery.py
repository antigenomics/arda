"""Build coordination must survive process death and serialize threads as well as processes."""
import multiprocessing as mp
import threading
import time

from arda import _locking


def _hold(lock, ready):
    with _locking.build_lock(lock, done=lambda: False):
        ready.send(True)
        time.sleep(30)


def test_killed_holder_releases_lock_without_manual_cleanup(tmp_path, monkeypatch):
    ctx = mp.get_context("spawn")
    reader, writer = ctx.Pipe(duplex=False)
    lock = tmp_path / ".build.lock"
    holder = ctx.Process(target=_hold, args=(lock, writer))
    holder.start()
    try:
        assert reader.poll(10), "holder did not start"
        assert reader.recv()
        holder.terminate()
        holder.join(10)
        assert not holder.is_alive()
        monkeypatch.setattr(_locking, "LOCK_TIMEOUT_S", 0.2)
        with _locking.build_lock(lock, done=lambda: False) as ours:
            assert ours
    finally:
        if holder.is_alive():
            holder.terminate()
            holder.join(10)


def test_threads_cannot_enter_the_same_build_simultaneously(tmp_path):
    lock = tmp_path / ".build.lock"
    second_started = threading.Event()
    second_inside = threading.Event()

    def second():
        second_started.set()
        with _locking.build_lock(lock, done=lambda: False):
            second_inside.set()

    with _locking.build_lock(lock, done=lambda: False):
        thread = threading.Thread(target=second)
        thread.start()
        assert second_started.wait(5)
        assert not second_inside.wait(0.05)
    thread.join(5)
    assert not thread.is_alive()
    assert second_inside.is_set()
