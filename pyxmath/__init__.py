class PT():
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __repr__(self):
        return repr([self.x, self.y])

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        else:
            return self.x == other[0] and self.y == other[1]

__all__ = ['PT']