from wait_quiet import wait_until_quiet


def run(loads, below=2.0, timeout_s=30):
    remaining = list(loads)
    slept = []
    quiet = wait_until_quiet(below, timeout_s, read_load=lambda: remaining.pop(0),
                             sleep=slept.append, interval_s=10)
    return quiet, len(slept)


def test_a_quiet_machine_does_not_wait():
    assert run([0.4]) == (True, 0)


def test_it_waits_until_the_load_falls():
    assert run([8.0, 3.1, 1.9]) == (True, 2)


def test_it_gives_up_after_the_timeout():
    assert run([8.0] * 10) == (False, 3)


def test_the_limit_itself_is_not_quiet():
    assert run([2.0, 1.99]) == (True, 1)
