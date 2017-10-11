class Test:
    def __init__(self):
        self.a = 1
        self.b = self.a
        
    def test(self):
        print 'self.a = ', self.a
        print 'self.b = ', self.b
        self.b = 20
        print 'self.a = ', self.a
        print 'self.b = ', self.b
        
T = Test()
T.test()
