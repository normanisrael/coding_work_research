from keras.applications import MobileNet
from keras import models
from keras import layers
from keras.preprocessing.image import ImageDataGenerator
from keras import optimizers

TRAIN_DIR = '/codes/data/train'
VAL_DIR = '/codes/data/validation'

conv_base = MobileNet(weights='imagenet',
                include_top=False,
                input_shape=(224, 224, 3))

model = models.Sequential()
model.add(conv_base)
model.add(layers.Flatten())
model.add(layers.Dropout(0.5))
model.add(layers.Dense(256, activation='relu'))
model.add(layers.Dense(1, activation='sigmoid'))

conv_base.trainable = False

train_datagen = ImageDataGenerator(
                                rescale=1./255,
                                #rotation_range=2.5,
                                width_shift_range=0.2,
                                height_shift_range=0.2,
                                shear_range=0.2,
                                zoom_range=0.2,
                                horizontal_flip=True,
                                vertical_flip=True,
                                fill_mode='nearest')

test_datagen = ImageDataGenerator(rescale=1./255)

train_generator = train_datagen.flow_from_directory(
                                                    TRAIN_DIR,
                                                    target_size=(224, 224),
                                                    batch_size=8,
                                                    class_mode='binary')

validation_generator = test_datagen.flow_from_directory(
                                                        VAL_DIR,
                                                        target_size=(224, 224),
                                                        batch_size=8,
                                                        class_mode='binary')

model.compile(loss='binary_crossentropy',
                optimizer=optimizers.RMSprop(lr=1e-5),
                metrics=['acc'])
                
history = model.fit_generator(
                            train_generator,
                            steps_per_epoch=100,
                            epochs=1000,
                            validation_data=validation_generator,
                            validation_steps=50)

conv_base.trainable = True
set_trainable = False
for layer in conv_base.layers:
    if layer.name == 'conv_pw_13':  #block5 initially (file 0)
        set_trainable = True
    if set_trainable:
        layer.trainable = True
    else:
        layer.trainable = False


model.compile(loss='binary_crossentropy',
                optimizer=optimizers.RMSprop(lr=1e-5),
                metrics=['acc'])

history = model.fit_generator(
                            train_generator,
                            steps_per_epoch=100,
                            epochs=1000,
                            validation_data=validation_generator,
                            validation_steps=50)
                    
training_accuracy = history.history['acc']
validation_accuracy = history.history['val_acc']
training_loss = history.history['loss']
validation_loss = history.history['val_loss']

epochs = range(1, len(training_accuracy) + 1)

model_diagnostics = open("TL_model_diagnostics_mobnet.dat", "a")
#model_diagnostics.write('{:<10s}{:10s}{:10s}{:15s}{:>15s} \n'.format('epoch', 'training accuracy', 'validation accuracy', 'training loss', 'validation loss'))
for i in range(0, len(epochs)):
    model_diagnostics.write('{:<3d}{:10f}{:10f}{:15.9f}{:>15.9f} \n'.format(epochs[i], training_accuracy[i], validation_accuracy[i], training_loss[i], validation_loss[i]))

model_diagnostics.close()
